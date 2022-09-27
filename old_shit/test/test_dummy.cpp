/**
 * Just a simple copy'n'paste template class to quickly add tests
 */

#include "Base.hh"

using namespace std;

void test_dummy() {
    ALEPH_ASSERT_EQUAL(42, 42)
}

int main(int argc, const char *argv[]) {
    test_dummy();
    return 0;
}