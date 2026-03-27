type MInt = i64;
type WolframLibraryData = *mut std::ffi::c_void;
type MArgument = *mut std::ffi::c_void;

#[no_mangle]
pub extern "C" fn WolframLibrary_getVersion() -> MInt {
    8
}

#[no_mangle]
pub extern "C" fn WolframLibrary_initialize(_lib_data: WolframLibraryData) -> i32 {
    0
}

#[no_mangle]
pub unsafe extern "C" fn lenia_simulate(
    _lib_data: WolframLibraryData,
    _argc: MInt,
    _args: *mut MArgument,
    _res: MArgument,
) -> i32 {
    0
}
