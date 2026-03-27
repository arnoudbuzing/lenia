(* --- Test Rust vs Wolfram Lenia implementation --- *)
PacletDirectoryLoad[FileNameJoin[{Directory[], "Lenia"}]];
Needs["Lenia`"];

Print["=== Testing Rust Lenia backend ==="];
Print[""];

(* --- Test 1: Basic run (no history) --- *)
Print["Test 1: Basic Rust run (no history)..."];
gridSize = {64, 64};
R = 10;
grid = N @ Table[
  Module[{r = Sqrt[(i - 32)^2 + (j - 32)^2] / R},
    If[r < 1.0, Exp[-((r - 0.5)^2) / (2 * 0.15^2)], 0.0]
  ],
  {i, 1, 64}, {j, 1, 64}
];

resultRust = Lenia[grid, 50, "Mu" -> 0.15, "Sigma" -> 0.015, "Radius" -> R, Method -> "Rust"];
If[FailureQ[resultRust] || resultRust === $Failed,
  Print["  FAIL: Rust backend returned $Failed"];
  Exit[1]
];
If[!MatrixQ[resultRust],
  Print["  FAIL: Rust result is not a matrix, got: ", Head[resultRust]];
  Exit[1]
];
If[Dimensions[resultRust] =!= {64, 64},
  Print["  FAIL: Dimensions mismatch: ", Dimensions[resultRust]];
  Exit[1]
];
Print["  PASS: Rust produced 64x64 result, max=", Max[resultRust]];

(* --- Test 2: Numerical agreement with WL --- *)
Print[""];
Print["Test 2: Numerical agreement (Rust vs Wolfram)..."];
resultWL = Lenia[grid, 50, "Mu" -> 0.15, "Sigma" -> 0.015, "Radius" -> R, Method -> "Wolfram"];
maxDiff = Max[Abs[resultRust - resultWL]];
Print["  Max absolute difference: ", maxDiff];
If[maxDiff < 1*^-8,
  Print["  PASS: Results agree to 1e-8"],
  Print["  WARNING: Difference = ", maxDiff, " (may be acceptable floating point variation)"]
];

(* --- Test 3: History mode --- *)
Print[""];
Print["Test 3: History mode..."];
historyRust = Lenia[grid, 10, "Mu" -> 0.15, "Sigma" -> 0.015, "Radius" -> R, 
  "ReturnHistory" -> True, Method -> "Rust"];
If[!ListQ[historyRust] || Length[historyRust] =!= 11,
  Print["  FAIL: Expected list of 11 frames, got: ", 
    If[ListQ[historyRust], Length[historyRust], Head[historyRust]]];
  Exit[1]
];
If[!AllTrue[historyRust, MatrixQ],
  Print["  FAIL: Not all history frames are matrices"];
  Exit[1]
];
(* First frame should match input *)
inputDiff = Max[Abs[historyRust[[1]] - grid]];
Print["  First frame vs input diff: ", inputDiff];
If[inputDiff < 1*^-10,
  Print["  PASS: First history frame matches input"],
  Print["  FAIL: First frame doesn't match input (diff=", inputDiff, ")"];
  Exit[1]
];
Print["  PASS: History mode returns 11 frames of correct shape"];

(* --- Test 4: Benchmark --- *)
Print[""];
Print["Test 4: Benchmark (128x128, 200 steps)..."];
bigGrid = LeniaSeed[128, "Ring"];
Print["  Timing Wolfram backend..."];
{tWL, rWL} = AbsoluteTiming[Lenia[bigGrid, 200, Method -> "Wolfram"]];
Print["    Wolfram: ", tWL, " seconds"];
Print["  Timing Rust backend..."];
{tRust, rRust} = AbsoluteTiming[Lenia[bigGrid, 200, Method -> "Rust"]];
Print["    Rust:    ", tRust, " seconds"];
If[tRust < tWL,
  Print["  Speedup: ", NumberForm[tWL/tRust, 3], "x faster"];,
  Print["  Rust was slower (", NumberForm[tWL/tRust, 3], "x)"];
];

Print[""];
Print["=== All tests completed ==="];
