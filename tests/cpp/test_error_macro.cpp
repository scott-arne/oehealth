// Simulates the ERROR macro that wingdi.h defines on Windows. It must come
// before any include so that this file is the first to include result.h.
#define ERROR 0
#include "oehealth/result.h"

#if !defined(ERROR) || ERROR != 0
#error "oehealth/result.h must restore the caller's ERROR macro"
#endif
#undef ERROR

#include <gtest/gtest.h>

namespace OEHEALTH::test {

TEST(ErrorMacro, ResultHeaderCompilesWithErrorDefined) {
    EXPECT_NE(Severity::ERROR, Severity::WARNING);
    EXPECT_NE(Status::ERROR, Status::FAIL);
}

} // namespace OEHEALTH::test
