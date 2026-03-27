(* Debug what FFI types are available *)
Print["$Version: ", $Version];
Print["$VersionNumber: ", $VersionNumber];

(* Try to understand the RawPointer syntax *)
Print[""];
Print["RawPointer with 2 args: ", RawPointer["Real64", 0]];

(* Try LibraryFunctionLoad instead - use the compilerDemoBase example *)
Print[""];
Print["--- Trying LibraryFunctionLoad approach ---"];

(* Check if our lib can be loaded via LibraryFunctionLoad with simple args *)
lib = "/Users/arnoudb/github/lenia/Lenia/LibraryResources/MacOSX-ARM64/liblenia_rs.dylib";

(* Try loading with NumericArray approach *)
Print[""];
Print["--- Trying NumericArray approach ---"];

(* Create a NumericArray and use RawMemoryExport *)
data = NumericArray[Table[0.0, 16], "Real64"];
Print["NumericArray: ", data];
ptr = RawMemoryExport[data];
Print["Exported ptr: ", ptr];
Print["Head: ", Head[ptr]];
