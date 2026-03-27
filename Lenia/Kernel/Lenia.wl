BeginPackage["Lenia`"];

Lenia::usage = "Lenia[grid, steps, opts] runs the Lenia continuous cellular automaton.";
LeniaStep::usage = "LeniaStep[grid, kernelFFT, mu, sigma, dt] performs one iteration using FFT convolution.";
LeniaKernel::usage = "LeniaKernel[gridSize, radius] creates a normalized Lenia kernel at the given grid size. LeniaKernel[gridSize, radius, type] uses a specific kernel type: \"Bump\", \"GaussianRing\", \"StepRing\", \"Polynomial\", \"SmoothLife\", or \"Sharp\".";
LeniaKernelFFT::usage = "LeniaKernelFFT[gridSize, radius] returns the FFT of the normalized Lenia kernel.";
LeniaGrowth::usage = "LeniaGrowth[n, mu, sigma] is the bell-shaped growth function.";
LeniaBlob::usage = "LeniaBlob[dims, center, radius] generates a localized density blob.";
LeniaSeed::usage = "LeniaSeed[n] generates an n\[Times]n grid with interesting initial conditions. LeniaSeed[n, type] uses a specific seed type: \"Ring\", \"MultiRing\", \"RandomOrganism\", \"Constellation\", or \"Asymmetric\".";

Lenia::norust = "Rust library not found at `1`. Falling back to Wolfram Language implementation.";
Lenia::rustfail = "Rust function call failed (code `1`). Falling back to Wolfram Language implementation.";

Begin["`Private`"];

(* --- Seed Generation: LeniaSeed --- *)

(* Helper: place a single ring organism on a grid *)
placeRing[grid_, center_, R_, peakR_ : 0.5, width_ : 0.15, amp_ : 1.0] := Module[{n1, n2, contribution},
  {n1, n2} = Dimensions[grid];
  contribution = N @ Table[
    Module[{r = Sqrt[(i - center[[1]])^2 + (j - center[[2]])^2] / R},
      If[r < 1.0, amp * Exp[-((r - peakR)^2) / (2 * width^2)], 0.0]
    ],
    {i, 1, n1}, {j, 1, n2}
  ];
  Clip[grid + contribution, {0.0, 1.0}]
];

(* Single ring centered on grid *)
leniaSeedRing[n_, R_] := placeRing[ConstantArray[0.0, {n, n}], {n/2, n/2}, R];

(* Two separate ring organisms *)
leniaSeedMultiRing[n_, R_] := Module[{grid, sep},
  grid = ConstantArray[0.0, {n, n}];
  sep = Round[R * 1.5];
  grid = placeRing[grid, {n/2, n/2 - sep}, R, 0.5, 0.15, 1.0];
  grid = placeRing[grid, {n/2, n/2 + sep}, R, 0.5, 0.15, 1.0];
  grid
];

(* A single ring with random perturbation for organic feel *)
leniaSeedRandomOrganism[n_, R_] := Module[{grid, noise},
  grid = placeRing[ConstantArray[0.0, {n, n}], {n/2, n/2}, R, 0.5, 0.18, 1.0];
  noise = RandomReal[{-0.15, 0.15}, {n, n}];
  Clip[grid + grid * noise, {0.0, 1.0}]
];

(* Three organisms in a triangular arrangement *)
leniaSeedConstellation[n_, R_] := Module[{grid, mid, sep},
  grid = ConstantArray[0.0, {n, n}];
  mid = n / 2;
  sep = Round[R * 2.0];
  grid = placeRing[grid, {mid - sep, mid}, R, 0.5, 0.15, 1.0];
  grid = placeRing[grid, {mid + Round[sep/2], mid - sep}, R, 0.5, 0.15, 0.9];
  grid = placeRing[grid, {mid + Round[sep/2], mid + sep}, R, 0.5, 0.15, 0.8];
  grid
];

(* Asymmetric shape: off-center ring with a bump *)
leniaSeedAsymmetric[n_, R_] := Module[{grid, offset},
  grid = ConstantArray[0.0, {n, n}];
  offset = Round[R * 0.3];
  grid = placeRing[grid, {n/2 - offset, n/2}, R, 0.5, 0.15, 1.0];
  grid = placeRing[grid, {n/2 + offset, n/2 + offset}, Round[R * 0.5], 0.4, 0.2, 0.7];
  grid
];

(* Public interface *)
LeniaSeed[n_Integer] := LeniaSeed[n, RandomChoice[{"Ring", "MultiRing", "RandomOrganism", "Constellation", "Asymmetric"}]];

LeniaSeed[n_Integer, "Ring"] := leniaSeedRing[n, Max[3, Round[n / 10]]];
LeniaSeed[n_Integer, "MultiRing"] := leniaSeedMultiRing[n, Max[3, Round[n / 10]]];
LeniaSeed[n_Integer, "RandomOrganism"] := leniaSeedRandomOrganism[n, Max[3, Round[n / 10]]];
LeniaSeed[n_Integer, "Constellation"] := leniaSeedConstellation[n, Max[3, Round[n / 10]]];
LeniaSeed[n_Integer, "Asymmetric"] := leniaSeedAsymmetric[n, Max[3, Round[n / 10]]];


(* --- Seed Generation --- *)
LeniaBlob[dims : {_Integer, _Integer}, center : {_?NumericQ, _?NumericQ}, radius_?NumericQ, amplitude : (_?NumericQ) : 1.0] :=
  N @ Table[
    amplitude * Exp[-((i - center[[1]])^2 + (j - center[[2]])^2) / radius^2],
    {i, 1, dims[[1]]}, {j, 1, dims[[2]]}
  ];

(* --- Growth Function --- *)
LeniaGrowth[n_, mu_, sigma_] := 2.0 * Exp[-((n - mu)^2) / (2.0 * sigma^2)] - 1.0;

(* --- Kernel Generation --- *)
(* Builds kernel at full grid size for FFT-based convolution.
   Following Bert Chan's reference: kernel is placed at grid center,
   then fftshifted so the center is at (0,0) for correct FFT convolution. *)

(* Core radial profiles: given normalized r in [0,1], return kernel value *)
leniaKernelCore["Bump", r_] := If[0.0 < r < 1.0, Exp[4.0 - 1.0 / (r * (1.0 - r))], 0.0];
leniaKernelCore["GaussianRing", r_] := If[0.0 < r < 1.0, Exp[-((r - 0.5)^2) / (2.0 * 0.15^2)], 0.0];
leniaKernelCore["StepRing", r_] := If[0.25 < r < 0.75, 1.0, 0.0];
leniaKernelCore["Polynomial", r_] := If[0.0 < r < 1.0, (4.0 * r * (1.0 - r))^4, 0.0];
leniaKernelCore["SmoothLife", r_] := If[0.0 < r < 1.0,
  0.5 * (1.0 + Tanh[20.0 * (r - 0.25)]) * 0.5 * (1.0 + Tanh[20.0 * (0.75 - r)]), 0.0];
leniaKernelCore["Sharp", r_] := If[0.0 < r < 1.0,
  Exp[-((r - 0.5)^2) / (2.0 * 0.05^2)], 0.0];

LeniaKernel[gridSize : {_Integer, _Integer}, radius_Integer, kernelType_String : "Bump"] :=
  Module[{mid, r, core, kernelNorm},
    mid = Ceiling[gridSize / 2];
    r = N @ Table[
      Sqrt[((i - mid[[1]])^2 + (j - mid[[2]])^2)] / radius,
      {i, 1, gridSize[[1]]}, {j, 1, gridSize[[2]]}
    ];
    core = Map[leniaKernelCore[kernelType, #] &, r, {2}];
    kernelNorm = core / Total[core, 2];
    RotateLeft[kernelNorm, mid - 1]
  ];

LeniaKernelFFT[gridSize : {_Integer, _Integer}, radius_Integer, kernelType_String : "Bump"] :=
  Fourier[LeniaKernel[gridSize, radius, kernelType], FourierParameters -> {1, -1}];

(* --- Single Step (FFT-based circular convolution) --- *)
LeniaStep[grid_?MatrixQ, kernelFFT_?MatrixQ, mu_, sigma_, dt_] := Module[{potential, change},
  potential = Re[InverseFourier[Fourier[grid, FourierParameters -> {1, -1}] * kernelFFT, FourierParameters -> {1, -1}]];
  change = LeniaGrowth[potential, mu, sigma];
  Clip[grid + dt * change, {0.0, 1.0}]
];

(* ================================================================== *)
(* --- Rust Backend via ForeignFunctionLoad --- *)
(* ================================================================== *)

$leniaRustLib = None;
$leniaRustFn = None;
$leniaRustFnHist = None;

$leniaPacletDir = DirectoryName[$InputFileName, 2];

leniaRustLibPath[] := Module[{ext},
  ext = Switch[$OperatingSystem, "MacOSX", ".dylib", "Unix", ".so", "Windows", ".dll"];
  FileNameJoin[{$leniaPacletDir, "LibraryResources", $SystemID, "liblenia_rs" <> ext}]
];

loadRustLibrary[] := Module[{libPath},
  libPath = leniaRustLibPath[];
  If[FileExistsQ[libPath],
    $leniaRustLib = libPath;
    True,
    Message[Lenia::norust, libPath];
    False
  ]
];

initRustFunction[] := Module[{},
  If[$leniaRustLib === None, loadRustLibrary[]];
  If[$leniaRustLib =!= None && $leniaRustFn === None,
    (* Non-history variant: returns NumericArray *)
    $leniaRustFn = LibraryFunctionLoad[$leniaRustLib, "lenia_simulate",
      {{Real, 2, "Constant"}, {Real, 2, "Constant"}, Integer, Real, Real, Real},
      NumericArray
    ];
    (* History variant: returns NumericArray *)
    $leniaRustFnHist = LibraryFunctionLoad[$leniaRustLib, "lenia_simulate_history",
      {{Real, 2, "Constant"}, {Real, 2, "Constant"}, Integer, Real, Real, Real},
      NumericArray
    ];
  ];
  $leniaRustFn =!= None
];

leniaRunRust[grid_?MatrixQ, kernel_?MatrixQ, steps_Integer, mu_, sigma_, dt_, returnHistory_] :=
  If[returnHistory,
    $leniaRustFnHist[N[grid], N[kernel], steps, N[mu], N[sigma], N[dt]]
  ,
    $leniaRustFn[N[grid], N[kernel], steps, N[mu], N[sigma], N[dt]]
  ];

(* ================================================================== *)
(* --- Main Loop --- *)
(* ================================================================== *)

Options[Lenia] = {
  "Mu" -> 0.15,
  "Sigma" -> 0.015,
  "DT" -> 0.1,
  "Radius" -> 13,
  "ReturnHistory" -> False,
  "LeniaKernel" -> "Bump",
  Method -> "Wolfram"
};

Lenia[initialGrid_?MatrixQ, steps_Integer, OptionsPattern[]] :=
  Module[{kernel, kernelFFT, grid, mu, sigma, dt, radius, returnHistory, method, kernelType},
    mu = OptionValue["Mu"];
    sigma = OptionValue["Sigma"];
    dt = OptionValue["DT"];
    radius = OptionValue["Radius"];
    returnHistory = TrueQ[OptionValue["ReturnHistory"]];
    method = OptionValue[Method];
    kernelType = OptionValue["LeniaKernel"];

    grid = N[initialGrid];

    (* Build kernel (shared by both backends) *)
    kernel = LeniaKernel[Dimensions[grid], radius, kernelType];

    (* Rust backend *)
    If[method === "Rust",
      If[initRustFunction[],
        Return[leniaRunRust[grid, kernel, steps, mu, sigma, dt, returnHistory]]
        (* else: Rust not available — fall through to Wolfram backend.
           Message already issued by initRustFunction. *)
      ]
    ];

    (* Wolfram Language backend *)
    kernelFFT = Fourier[kernel, FourierParameters -> {1, -1}];

    If[returnHistory,
      NumericArray[NestList[LeniaStep[#, kernelFFT, mu, sigma, dt] &, grid, steps], "Real64"]
    ,
      NumericArray[Nest[LeniaStep[#, kernelFFT, mu, sigma, dt] &, grid, steps], "Real64"]
    ]
  ];

End[];
EndPackage[];
