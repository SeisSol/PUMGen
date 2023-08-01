#ifndef APF_CONVERT_WRAPPER_H
#define APF_CONVERT_WRAPPER_H

#include <apfConvert.h>
#include <type_traits>

// This is a small fix, because some time after PUMI (i.e. SCOREC/core) v2.2.7,
// the element ID was switched from a mere int to a named type apf::GiD, a typedef to long.
// Thus, PUMgen did not compile anymore, but to my knowledge, there was no flag to check the version.
// (of course, I could have missed some)
// To make PUMgen work with newer versions of PUMI, we therefore read out the type of the only function
// that actually was affected by the switch to apf::GiD and use that one in PUMgen.

// inspired by https://stackoverflow.com/a/70954691

namespace internal {
template<typename T>
struct GetSecondFunctionArgument {};

template<typename F, typename Arg1, typename Arg2, typename ...ArgRest>
struct GetSecondFunctionArgument<F(Arg1, Arg2, ArgRest...)> {
    using Type = Arg2;
};

using ElementIDConstPtr = GetSecondFunctionArgument<decltype(apf::construct)>::Type;
using ElementID = std::remove_const_t<std::remove_pointer_t<ElementIDConstPtr>>;

}

using ElementID = internal::ElementID;

#endif
