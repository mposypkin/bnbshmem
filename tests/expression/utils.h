#ifndef UTILS__HPP
#define UTILS__HPP

#include <atomic>

namespace snowgoose {
namespace expression {

class Utils {
public :
    static int getUid();
};

int Utils::getUid() {
    static std::atomic<std::uint32_t> uid { 0 };
    return ++uid;
}

}}

#endif
