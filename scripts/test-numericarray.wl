PacletDirectoryLoad[FileNameJoin[{Directory[], "Lenia"}]];
Needs["Lenia`"];

Print["=== Testing NumericArray Return Type ==="];

grid = LeniaSeed[32, "Ring"];

(* Test 1: Rust backend (no history) *)
Print["Test 1: Rust backend (no history)..."];
res1 = Lenia[grid, 10, Method -> "Rust"];
If[Head[res1] === NumericArray,
  Print["  PASS: Head is NumericArray. Dimensions: ", Dimensions[res1]],
  Print["  FAIL: Head is ", Head[res1]]
];

(* Test 2: Rust backend (history) *)
Print["Test 2: Rust backend (history)..."];
res2 = Lenia[grid, 5, "ReturnHistory" -> True, Method -> "Rust"];
If[Head[res2] === NumericArray,
  Print["  PASS: Head is NumericArray. Dimensions: ", Dimensions[res2]],
  Print["  FAIL: Head is ", Head[res2]]
];

(* Test 3: Wolfram backend (no history) *)
Print["Test 3: Wolfram backend (no history)..."];
res3 = Lenia[grid, 10, Method -> "Wolfram"];
If[Head[res3] === NumericArray,
  Print["  PASS: Head is NumericArray. Dimensions: ", Dimensions[res3]],
  Print["  FAIL: Head is ", Head[res3]]
];

(* Test 4: Wolfram backend (history) *)
Print["Test 4: Wolfram backend (history)..."];
res4 = Lenia[grid, 5, "ReturnHistory" -> True, Method -> "Wolfram"];
If[Head[res4] === NumericArray,
  Print["  PASS: Head is NumericArray. Dimensions: ", Dimensions[res4]],
  Print["  FAIL: Head is ", Head[res4]]
];

Print["Done."];
