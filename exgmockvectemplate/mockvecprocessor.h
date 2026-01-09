// MockVecProcessor.hpp
#pragma once
#include <gmock/gmock.h>
#include "vecwrapper.h"

template <typename T>
class IVecProcessor {
public:
    virtual ~IVecProcessor() = default;
    virtual void process(VecWrapper<T>& v) = 0;
};

template <typename T>
class MockVecProcessor : public IVecProcessor<T> {
public:
    MOCK_METHOD(void, process, (VecWrapper<T>& v), (override));
};