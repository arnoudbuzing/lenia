#include <stdio.h>
#include <stddef.h>
#include "/Applications/Wolfram/15.0/Wolfram.app/Contents/SystemFiles/IncludeFiles/C/WolframLibrary.h"

int main() {
    struct st_WolframLibraryData lib;
    printf("MTENSOR_NEW_OFFSET: %lu\n", (unsigned long)((char*)&lib.MTensor_new - (char*)&lib) / sizeof(void*));
    printf("MTENSOR_GETRANK_OFFSET: %lu\n", (unsigned long)((char*)&lib.MTensor_getRank - (char*)&lib) / sizeof(void*));
    printf("MTENSOR_GETDIMENSIONS_OFFSET: %lu\n", (unsigned long)((char*)&lib.MTensor_getDimensions - (char*)&lib) / sizeof(void*));
    printf("MTENSOR_GETREALDATA_OFFSET: %lu\n", (unsigned long)((char*)&lib.MTensor_getRealData - (char*)&lib) / sizeof(void*));
    printf("MTENSOR_GETCOMPLEXDATA_OFFSET: %lu\n", (unsigned long)((char*)&lib.MTensor_getComplexData - (char*)&lib) / sizeof(void*));
    printf("VERSION_NUMBER_OFFSET: %lu\n", (unsigned long)((char*)&lib.VersionNumber - (char*)&lib) / sizeof(void*));
    return 0;
}
