/* As far as I know, there is no proper cname mechanism for this, so
 * we have to manually unmangle this Cython name. */
#define ObjectPoolObject __pyx_obj_4sage_3ext_4pool_ObjectPool


static long pools_used = 0;


#define DEFINE_POOL_NUMBER(i) \
static struct ObjectPoolObject* pool_##i; \
static PyObject* pool_##i##_new(PyTypeObject* cls, PyObject* args, PyObject* kwds) {return ObjectPool_new(pool_##i, (PyObject*)cls, args, kwds);} \
static void pool_##i##_dealloc(PyObject* obj) {ObjectPool_dealloc(pool_##i, obj);}

/* The number of DEFINE_POOL_NUMBER() invocations must equal the number of
 # CASE() invocations below. */
DEFINE_POOL_NUMBER(0);
DEFINE_POOL_NUMBER(1);
DEFINE_POOL_NUMBER(2);
DEFINE_POOL_NUMBER(3);


static CYTHON_INLINE void set_functions(PyTypeObject* t, newfunc tp_new, destructor tp_dealloc)
{
    t->tp_new = tp_new;
    t->tp_dealloc = tp_dealloc;
}


#define CASE(i) \
    case i: pool_##i = p; set_functions(t, pool_##i##_new, pool_##i##_dealloc); return pools_used = i+1;

static long _attach_pool(PyTypeObject* t, struct ObjectPoolObject* p)
{
    /* Put the pool in the type dict as the __objectpool attribute.
     * This also ensures that the pool won't be garbage collected. */
    if (PyDict_SetItemString(t->tp_dict, "__ObjectPool", (PyObject*)p) == -1) return -1;

    switch(pools_used)
    {
        CASE(0);
        CASE(1);
        CASE(2);
        CASE(3);
    }

    PyErr_Format(PyExc_OverflowError, "at most %ld object pools are supported per Cython module", pools_used);
    return -1;
}


static CYTHON_INLINE void Py_DECREF_KEEP(PyObject* obj)
{
    --Py_REFCNT(obj);
}
