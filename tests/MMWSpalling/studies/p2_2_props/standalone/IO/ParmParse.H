// Minimal stand-in for IO/ParmParse.H. Grimsel::Parse is never CALLED offline
// (the dump default-constructs), but it must COMPILE, so these signatures
// mirror the real ones and do nothing.
#ifndef IO_PARMPARSE_H
#define IO_PARMPARSE_H
#include <string>
namespace IO {
class ParmParse
{
public:
    ParmParse() = default;
    explicit ParmParse(std::string) {}
    // ⚠ ONE overload only. A second taking (std::string, T&, double) is
    // ambiguous against this one for every T = double call, which is every
    // call Grimsel::Parse makes.
    template <class T> void query_default(std::string, T &, T) {}
    template <class C> void queryclass(std::string, C &) {}
};
}
#endif
