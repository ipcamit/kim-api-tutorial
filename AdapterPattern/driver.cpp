#include "driver.h"

// model driver
class MyDriver {
public:
    MyDriver(int p) : p_(p) {}

    /* NON-static member — note implicit `this`               */
    int compute(int x) { return this->p_ + x; } // this = self in Python

    /* STATIC member — *no* implicit `this`, plain C signature */
    static int compute_adapter(void* ctx, int x) // equivalent to Compute functions in driver
    {
        auto* self = static_cast<MyDriver*>(ctx);  /* recover the object */
        return self->compute(x);                   /* delegate */
    }

private:
    int p_;
};

/* ----------  C-visible function  ---------- */
compute_fn make_driver(int parameter, void** out_ctx)
{
    MyDriver* obj = new MyDriver(parameter);   /* heap-allocate object  */
    *out_ctx = obj;                            /* give C code the handle */
    return &MyDriver::compute_adapter;         /* return C-style fn ptr */
}
