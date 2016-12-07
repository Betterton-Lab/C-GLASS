#ifndef _SIMCORE_TEST_XLINK_KMC_MOVES_H_
#define _SIMCORE_TEST_XLINK_KMC_MOVES_H_

#include "xlink_kmc.h"
#include "space.h"
#include "test_module_base.h"

#include <functional>

class BrRod;
class Xlink;

class TestXlinkKMCMoves : public TestModuleBase, public XlinkKMC {
  public:
    TestXlinkKMCMoves() : TestModuleBase(), XlinkKMC() {}
    virtual ~TestXlinkKMCMoves() {}

    virtual void InitTestModule(const std::string& filename);
    virtual void RunTests();

  protected:

    bool TestKMC_0_1(int test_num);
    bool TestKMC_1_0(int test_num);
    bool TestKMC_1_2(int test_num);

    system_parameters params_sub_;
    SpaceProperties space_sub_;
};

#endif