#include "oehealth/result.h"

namespace OEHEALTH {

Status worst_status(Status a, Status b) {
    return static_cast<int>(a) >= static_cast<int>(b) ? a : b;
}

} // namespace OEHEALTH
