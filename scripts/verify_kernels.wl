PacletDirectoryLoad[FileNameJoin[{Directory[], "Lenia"}]];
Needs["Lenia`"];

grid = LeniaSeed[64, "Ring"];
kernels = {"Bump", "GaussianRing", "StepRing", "Polynomial", "SmoothLife", "Sharp"};

Do[
  Print["Testing kernel: ", k];
  res = Lenia[grid, 10, "LeniaKernel" -> k, Method -> "Rust"];
  If[MatrixQ[res], 
    Print["  PASS: Max=", Max[res]],
    Print["  FAIL: Got ", Head[res]]
  ],
  {k, kernels}
];
Print["All kernel types verified."];
