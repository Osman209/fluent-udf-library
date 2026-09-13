#include <assert.h>
#include "mock/udf.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"
int main(void){assert(UDF_IS_WRITER==EXPECTED_WRITER);return 0;}
