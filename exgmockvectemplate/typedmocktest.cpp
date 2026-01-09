// TypedMockTest.cpp
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "vecwrapper.h"
#include "mockvecprocessor.h"

using ::testing::Invoke;
using ::testing::_;

// Define the types to test
using MyTypes = ::testing::Types<int, double>;
template <typename T>
class VecWrapperTypedTest : public ::testing::Test {};
TYPED_TEST_SUITE(VecWrapperTypedTest, MyTypes);

// The actual typed test
TYPED_TEST(VecWrapperTypedTest, MockedProcessingAddsOneElement)
{
    using T = TypeParam;

    VecWrapper<T> vec;
    vec.push(T{});

    MockVecProcessor<T> mock;

    // Mock behaviour: push one extra element during processing
    EXPECT_CALL(mock, process(_))
        .WillOnce(Invoke([](VecWrapper<T>& v){
            v.push(T{});
        }));

    mock.process(vec);

    EXPECT_EQ(vec.size(), 2u);
}