#pragma once

#ifdef __CUDACC__
    #define CUDA_FUNC __device__ __host__
#else
    #define CUDA_FUNC
#endif


CUDA_FUNC inline float3 make_float3(const float4 &v) {
	return make_float3(v.x, v.y, v.z);
}

CUDA_FUNC inline float4 make_float4(const float3 &v, float w = 0.0f) {
	return make_float4(v.x, v.y, v.z, w);
}

// Math overloading

CUDA_FUNC inline float4 operator+(const float4& a, const float4& b) {
	return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

CUDA_FUNC inline float3 operator+(const float3& a, const float3& b) {
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

CUDA_FUNC inline float4 operator-(const float4& a, const float4& b) {
	return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

CUDA_FUNC inline float3 operator-(const float3& a, const float3& b) {
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

CUDA_FUNC inline void operator+=(float4& a, const float4& b) {
	a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

CUDA_FUNC inline void operator+=(float3& a, const float3& b) {
	a.x += b.x; a.y += b.y; a.z += b.z;
}

CUDA_FUNC inline void operator-=(float4& a, const float4& b) {
	a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}

CUDA_FUNC inline void operator-=(float3& a, const float3& b) {
	a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

CUDA_FUNC inline void operator/=(float4& a, float b) {
	a.x /= b; a.y /= b; a.z /= b; a.w /= b;
}

CUDA_FUNC inline void operator/=(float3& a, float b) {
	a.x /= b; a.y /= b; a.z /= b;
}

CUDA_FUNC inline float4 operator/(const float4& a, float b) {
	return make_float4(a.x / b, a.y / b, a.z / b, a.w / b);
}

CUDA_FUNC inline float3 operator/(const float3& a, float b) {
	return make_float3(a.x / b, a.y / b, a.z / b);
}

CUDA_FUNC inline void operator*=(float4& a, float b) {
	a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

CUDA_FUNC inline void operator*=(float3& a, float b) {
	a.x *= b; a.y *= b; a.z *= b;
}

CUDA_FUNC inline float4 operator*(const float4& a, float b) {
	return make_float4(a.x * b, a.y * b, a.z * b, a.w * b);
}

CUDA_FUNC inline float4 operator*(float b, const float4& a) {
	return a * b;
}

CUDA_FUNC inline float3 operator*(const float3& a, float b) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

CUDA_FUNC inline float3 operator*(float b, const float3& a) {
	return a * b;
}

// 3D Vector operations

template <typename T>
CUDA_FUNC inline float dot(const T &a, const T &b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

CUDA_FUNC inline float4 cross(const float4 &a, const float4 &b) {
	return make_float4(a.y*b.z-b.y*a.z, a.z*b.x-b.z*a.x, a.x*b.y-b.x*a.y, 0.0f);
}

CUDA_FUNC inline float3 cross(const float3 &a, const float3 &b) {
	return make_float3(a.y*b.z-b.y*a.z, a.z*b.x-b.z*a.x, a.x*b.y-b.x*a.y);
}

template <typename T>
CUDA_FUNC inline float abs2(const T &a) {
	return a.x * a.x + a.y * a.y + a.z * a.z;
}

CUDA_FUNC inline float abs(const float4 &a) {
	return sqrtf(abs2(a));
}

CUDA_FUNC inline float abs(const float3 &a) {
	return sqrtf(abs2(a));
}

template <typename T>
CUDA_FUNC inline float normalize(T &a) {
	float v = abs(a);
	a /= v;
	return v;
}

template <typename T>
CUDA_FUNC inline float normalize_s(T &a) { // Safe normalize
	float v = abs(a);
	if (v < 1e-5) // Chosen arbitrary. It just seemed suitable
		a = make_float4(1.0f, 0.0f, 0.0f, 0.0f); // Chosen arbitrary
	else
		a /= v;
	return v;
}

// Aargh, `isnan' is actually a macro, so no overloading :-(
// I've read somewhere that it will be made template in C++0x, but who cares
// Altogether, this piece of code sucks, but I'm lazy
CUDA_FUNC inline int my_isnan(const float4 &a) {
	return (isnan(a.x) || isnan(a.y) || isnan(a.z) || isnan(a.w));
}

CUDA_FUNC inline int my_isnan(const float3 &a) {
	return (isnan(a.x) || isnan(a.y) || isnan(a.z));
}

CUDA_FUNC inline int my_isnan(const float2 &a) {
		return (isnan(a.x) || isnan(a.y));
}

CUDA_FUNC inline int my_isnan(const float &a) {
	return isnan(a);
}

inline std::ostream& operator << (std::ostream& os, int4 v) {
	return os << "[" << v.x << ", " << v.y << ", " << v.z << "; " << v.w << "]";
}

inline std::ostream& operator << (std::ostream& os, float3 v) {
	return os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
}

inline std::ostream& operator << (std::ostream& os, float4 v) {
	return os << "[" << v.x << ", " << v.y << ", " << v.z << "; " << v.w << "]";
}

