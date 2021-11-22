//
// Created by kulaginasv on 22.11.21.
//

#include "../googletest-main/googletest/include/gtest/gtest.h"
#include "../include/tree.h"

namespace {

    class TreeTest : public ::testing::Test {
    protected:

        TreeTest() {
            // You can do set-up work for each test here.
        }

        ~TreeTest() override {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        void SetUp() override {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        void TearDown() override {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

    };

// Tests that the Foo::Bar() method does Abc.
    TEST_F(TreeTest, MethodBarDoesAbc) {
    Tree *f;
    EXPECT_EQ(f->getSize(), 0);
}

// Tests that Foo does Xyz.
TEST_F(TreeTest, DoesXyz) {
// Exercises the Xyz feature of Foo.
}

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
