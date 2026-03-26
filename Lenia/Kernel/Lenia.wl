BeginPackage["Lenia`"];

Lenia::usage = "Lenia[grid, steps, opts] runs the Lenia continuous cellular automaton.";
LeniaStep::usage = "LeniaStep[grid, kernelFFT, mu, sigma, dt] performs one iteration using FFT convolution.";
LeniaKernel::usage = "LeniaKernel[gridSize, radius] creates a normalized Lenia kernel at the given grid size.";
LeniaKernelFFT::usage = "LeniaKernelFFT[gridSize, radius] returns the FFT of the normalized Lenia kernel.";
LeniaGrowth::usage = "LeniaGrowth[n, mu, sigma] is the bell-shaped growth function.";
LeniaBlob::usage = "LeniaBlob[dims, center, radius] generates a localized density blob.";
LeniaSeed::usage = "LeniaSeed[n] generates an n\[Times]n grid with interesting initial conditions. LeniaSeed[n, type] uses a specific seed type: \"Ring\", \"MultiRing\", \"RandomOrganism\", \"Constellation\", or \"Asymmetric\".";

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

(* Two concentric rings *)
leniaSeedMultiRing[n_, R_] := Module[{grid},
  grid = ConstantArray[0.0, {n, n}];
  grid = placeRing[grid, {n/2, n/2}, R, 0.35, 0.1, 0.8];
  grid = placeRing[grid, {n/2, n/2}, R, 0.7, 0.08, 0.6];
  grid
];

(* A single ring with random perturbation for organic feel *)
leniaSeedRandomOrganism[n_, R_] := Module[{grid, noise},
  grid = placeRing[ConstantArray[0.0, {n, n}], {n/2, n/2}, R, 0.5, 0.18, 1.0];
  noise = RandomReal[{-0.15, 0.15}, {n, n}];
  Clip[grid + grid * noise, {0.0, 1.0}]
];

(* Multiple small organisms scattered across the grid *)
leniaSeedConstellation[n_, R_] := Module[{grid, nOrganisms, centers, smallR},
  grid = ConstantArray[0.0, {n, n}];
  nOrganisms = RandomInteger[{3, 6}];
  smallR = Max[3, Round[R * 0.6]];
  centers = Table[
    {RandomInteger[{Round[n * 0.2], Round[n * 0.8]}],
     RandomInteger[{Round[n * 0.2], Round[n * 0.8]}]},
    {nOrganisms}
  ];
  Do[
    grid = placeRing[grid, c, smallR, RandomReal[{0.3, 0.6}], RandomReal[{0.1, 0.2}], RandomReal[{0.5, 1.0}]],
    {c, centers}
  ];
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
LeniaKernel[gridSize : {_Integer, _Integer}, radius_Integer] := Module[{mid, r, core, kernelNorm},
  mid = Ceiling[gridSize / 2];
  (* Normalized distance: r = distance_from_center / R *)
  r = N @ Table[
    Sqrt[((i - mid[[1]])^2 + (j - mid[[2]])^2)] / radius,
    {i, 1, gridSize[[1]]}, {j, 1, gridSize[[2]]}
  ];
  (* Exponential bump kernel: nonzero only for 0 < r < 1 *)
  core = Map[If[0.0 < # < 1.0, Exp[4.0 - 1.0 / (# * (1.0 - #))], 0.0] &, r, {2}];
  (* Normalize to sum to 1 *)
  kernelNorm = core / Total[core, 2];
  (* Shift center to (0,0) corner for FFT convolution *)
  RotateLeft[kernelNorm, mid - 1]
];

LeniaKernelFFT[gridSize : {_Integer, _Integer}, radius_Integer] :=
  Fourier[LeniaKernel[gridSize, radius], FourierParameters -> {1, -1}];

(* --- Single Step (FFT-based circular convolution) --- *)
LeniaStep[grid_?MatrixQ, kernelFFT_?MatrixQ, mu_, sigma_, dt_] := Module[{potential, change},
  potential = Re[InverseFourier[Fourier[grid, FourierParameters -> {1, -1}] * kernelFFT, FourierParameters -> {1, -1}]];
  change = LeniaGrowth[potential, mu, sigma];
  Clip[grid + dt * change, {0.0, 1.0}]
];

(* --- Main Loop --- *)
Options[Lenia] = {
  "Mu" -> 0.15,
  "Sigma" -> 0.015,
  "DT" -> 0.1,
  "Radius" -> 13,
  "ReturnHistory" -> False
};

Lenia[initialGrid_?MatrixQ, steps_Integer, OptionsPattern[]] := Module[{kernelFFT, grid, mu, sigma, dt, radius},
  mu = OptionValue["Mu"];
  sigma = OptionValue["Sigma"];
  dt = OptionValue["DT"];
  radius = OptionValue["Radius"];
  
  grid = N[initialGrid];
  kernelFFT = LeniaKernelFFT[Dimensions[grid], radius];
  
  If[OptionValue["ReturnHistory"],
    NestList[LeniaStep[#, kernelFFT, mu, sigma, dt] &, grid, steps]
  ,
    Nest[LeniaStep[#, kernelFFT, mu, sigma, dt] &, grid, steps]
  ]
];

End[];
EndPackage[];
