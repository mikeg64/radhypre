
// VecWrapper.hpp
#pragma once
#include <vector>
#include <stdexcept>

template <typename T>
class VecWrapper {
public:
    VecWrapper() = default;

    explicit VecWrapper(std::initializer_list<T> init)
        : data_(init) {}

    void push(const T& value) {
        data_.push_back(value);
    }

    std::size_t size() const {
        return data_.size();
    }

    const T& operator[](std::size_t i) const {
        return data_[i];
    }

    T& operator[](std::size_t i) {
        return data_[i];
    }

    // Vector addition
    VecWrapper<T> operator+(const VecWrapper<T>& other) const {
        if (data_.size() != other.data_.size()) {
            throw std::runtime_error("Vector sizes do not match");
        }

        VecWrapper<T> result;
        result.data_.resize(data_.size());

        for (std::size_t i = 0; i < data_.size(); i++) {
            result.data_[i] = data_[i] + other.data_[i];
        }
        return result;
    }

private:
    std::vector<T> data_;
};



