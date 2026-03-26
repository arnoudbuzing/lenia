(* --- Lenia Unit Tests --- *)

Needs["Lenia`"];

VerificationTest[
  Total[LeniaKernel[{32, 32}, 5], 2],
  1.0,
  {SameTest -> (Abs[#1 - #2] < 10^-10 &)},
  TestID -> "KernelNormalization"
]

VerificationTest[
  LeniaGrowth[0.15, 0.15, 0.015],
  1.0,
  TestID -> "GrowthAtMu"
]

VerificationTest[
  Quiet[LeniaGrowth[1.0, 0.15, 0.015]],
  -1.0,
  {SameTest -> (Abs[#1 - #2] < 10^-10 &)},
  TestID -> "GrowthFarFromMu"
]

(* Verify FFT convolution correctness: constant field -> constant potential *)
VerificationTest[
  Module[{gridSize = {32, 32}, kFFT, potential},
    kFFT = LeniaKernelFFT[gridSize, 5];
    potential = Re[InverseFourier[Fourier[ConstantArray[0.5, gridSize], FourierParameters -> {1, -1}] * kFFT, FourierParameters -> {1, -1}]];
    potential[[16, 16]]
  ],
  0.5,
  {SameTest -> (Abs[#1 - #2] < 10^-6 &)},
  TestID -> "FFTConvolutionConstant"
]

(* LeniaStep preserves grid dimensions *)
VerificationTest[
  Module[{gridSize = {32, 32}, grid, kFFT},
    grid = LeniaBlob[gridSize, {16, 16}, 5, 0.5];
    kFFT = LeniaKernelFFT[gridSize, 5];
    Dimensions[LeniaStep[grid, kFFT, 0.15, 0.015, 0.1]]
  ],
  {32, 32},
  TestID -> "StepShapePreservation"
]

(* Lenia main function runs without error *)
VerificationTest[
  MatrixQ[Lenia[LeniaBlob[{32, 32}, {16, 16}, 5, 0.5], 10, "Radius" -> 5]],
  True,
  TestID -> "LeniaBasicRun"
]

(* Simulation does NOT decay to zero with ring-shaped seed *)
VerificationTest[
  Module[{gridSize = {64, 64}, R = 10, grid, result},
    grid = N @ Table[
      Module[{r = Sqrt[(i - 32)^2 + (j - 32)^2] / R},
        If[r < 1.0, Exp[-((r - 0.5)^2) / (2 * 0.15^2)], 0.0]
      ],
      {i, 1, 64}, {j, 1, 64}
    ];
    result = Lenia[grid, 100, "Mu" -> 0.15, "Sigma" -> 0.015, "Radius" -> R];
    Max[result] > 0.01
  ],
  True,
  TestID -> "SimulationSurvival"
]

(* LeniaBlob shape and amplitude *)
VerificationTest[
  Dimensions[LeniaBlob[{50, 50}, {25, 25}, 5]],
  {50, 50},
  TestID -> "BlobShape"
]

VerificationTest[
  Max[LeniaBlob[{50, 50}, {25, 25}, 5, 1.0]],
  1.0,
  TestID -> "BlobAmplitude"
]
