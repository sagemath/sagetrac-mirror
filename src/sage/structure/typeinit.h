/*****************************************************************************
*       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 2 of the License, or
* (at your option) any later version.
*                  http://www.gnu.org/licenses/
*****************************************************************************/

/*
 * This function calls PyType_Ready(t) and then calls t.__typeinit__(t)
 * as if that was a class method. The __typeinit__ method can then be
 * used for example to make changes to the tp_... slots.
 */
static CYTHON_INLINE int Sage_PyType_Ready(PyTypeObject* t)
{
    int r = PyType_Ready(t);
    if (r < 0) return r;

    PyObject* typeinit;
    typeinit = PyObject_GetAttrString((PyObject*)t, "__typeinit__");
    if (!typeinit)
    {
        /* No __typeinit__ method: no problem, go on... */
        PyErr_Clear();
        return 0;
    }

    /* Now typeinit is a method, but we want to call it as if it was
     * a classmethod. Some hackery is needed here... */
    PyMethodDef* meth = ((PyMethodDescrObject*)typeinit)->d_method;
    Py_DECREF(typeinit);

    /* This must be a METH_NOARGS method */
    if (meth->ml_flags != METH_NOARGS)
    {
        PyErr_SetString(PyExc_TypeError,
                "__typeinit__ must be an ordinary method without arguments");
    }

    PyObject* ret = meth->ml_meth((PyObject*)t, NULL);
    if (!ret) return -1;

    Py_DECREF(ret);
    return 0;
}

/* Use the above function in Cython code instead of the default
 * PyType_Ready() function */
#define PyType_Ready(t)  Sage_PyType_Ready(t)
