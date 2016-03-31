struct __stack_element {
	int id;
	double tmin, tmax;
};

class __stack {
public:
	__stack_element stack[32];
	int count;
	__device__ __stack();
	__device__ void push(int id, double tmin, double tmax);
	__device__ __stack_element pop();
	__device__ bool empty();
};

__device__ __stack::__stack() : count(0) {}

__device__ void __stack::push(int id, double tmin, double tmax) {
	this->stack[count].id = id;
	this->stack[count].tmin = tmin;
	this->stack[count].tmax = tmax;
	count++;
}

__device__ __stack_element __stack::pop() {
	count--;
	__stack_element se;
	se.id = this->stack[count].id;
	se.tmin = this->stack[count].tmin;
	se.tmax = this->stack[count].tmax;
	return se;
}

__device__ bool __stack::empty() {
	return this->count == 0;
}
