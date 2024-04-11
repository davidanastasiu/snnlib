/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include <cstdio>

namespace snnlib {

    /* Macros needed by the qsort routines */
    /* Discontinue quicksort algorithm when partition gets below this size.
       This particular magic number was chosen to work best on a Sun 4/260. */
    #define DA_QSORT_MAX_THRESH 8
    /* The next 4 #defines implement a very fast in-line stack abstraction. */
    #define DA_QSORT_STACK_SIZE     (8 * sizeof(size_t))
    /* Swap two items pointed to by A and B using temporary buffer t. */
    // #define DA_QSORT_SWAP(a, b, t) ((void)((t = *a), (*a = *b), (*b = t)))
    // #define DA_QSORT_PUSH(top, low, high) (((top->_lo = (low)), (top->_hi = (high)), ++top))
    // #define DA_QSORT_POP(low, high, top)    ((--top, (low = top->_lo), (high = top->_hi)))
    /* qsort routine comparison for ascending direction */
    #define DA_QSORT_LT(a, b) ((*a) < (*b))
    /* qsort routine comparison decreasing direction */
    #define DA_QSORT_GT(a, b) ((*a) > (*b))



    /**
     * Sort key-value pairs in the keys and values arrays in ascending oder by value.
     * In-line qsort implementation.    Differs from traditional qsort() routine in that
     * it does not call a function to swap elements. Moreover, it sorts two arrays
     * simultaneously based on the values in the vals array and swaps the associated
     * key with each swapped value.
     *
     * @param n     Size of keys and vals arrays.
     * @param keys    Array of keys.
     * @param vals    Array of values.
     */
    /* Copyright (C) 1991, 1992, 1996, 1997, 1999 Free Software Foundation, Inc.
        This file is part of the GNU C Library.
        Written by Douglas C. Schmidt (schmidt@ics.uci.edu).

        The GNU C Library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Lesser General Public
        License as published by the Free Software Foundation; either
        version 2.1 of the License, or (at your option) any later version.

        The GNU C Library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the GNU
        Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public
        License along with the GNU C Library; if not, write to the Free
        Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
        02111-1307 USA.    */
    template<class K, class V>
    void kvsorti(const size_t n, K *const keys, V *const vals) {
        {
            V *const _base = (vals);
            const size_t _elems = (n);
            K ktmp;
            V vtmp;
            if (_elems == 0) {
                return;
            }
            if (_elems > DA_QSORT_MAX_THRESH) {
                V *_lo = _base;
                V *_hi = _lo + _elems - 1;
                struct {
                    V *_hi;
                    V *_lo;
                } _stack[DA_QSORT_STACK_SIZE], *_top = _stack + 1;
                while ((_stack < _top)) {
                    V *_left_ptr;
                    V *_right_ptr;
                    /* Select median value from among LO, MID, and HI. Rearrange    \
                         LO and HI so the three values are sorted. This lowers the    \
                         probability of picking a pathological pivot value and        \
                         skips a comparison for both the LEFT_PTR and RIGHT_PTR in    \
                         the while loops. */
                    V *_mid = _lo + ((_hi - _lo) >> 1);
                    if (((*_mid) < (*_lo))) {    /* if (QSORT_LT (_mid, _lo)) */
                        /* QSORT_SWAP (_mid, _lo, _hold); */
                        vtmp = *_mid;
                        ktmp = keys[_mid - vals];
                        *_mid = *_lo;
                        keys[_mid - vals] = keys[_lo - vals];
                        *_lo = vtmp;
                        keys[_lo - vals] = ktmp;
                        /* end QSORT_SWAP */
                    }
                    if (((*_hi) < (*_mid))) {    /* if (QSORT_LT (_hi, _mid)) */
                        /* QSORT_SWAP (_mid, _hi, _hold); */
                        vtmp = *_mid;
                        ktmp = keys[_mid - vals];
                        *_mid = *_hi;
                        keys[_mid - vals] = keys[_hi - vals];
                        *_hi = vtmp;
                        keys[_hi - vals] = ktmp;
                        /* end QSORT_SWAP */
                    } else {
                        goto _jump_over;
                    }
                    if (((*_mid) < (*_lo))) {    /* if (QSORT_LT (_mid, _lo)) */
                        /* QSORT_SWAP (_mid, _lo, _hold); */
                        vtmp = *_mid;
                        ktmp = keys[_mid - vals];
                        *_mid = *_lo;
                        keys[_mid - vals] = keys[_lo - vals];
                        *_lo = vtmp;
                        keys[_lo - vals] = ktmp;
                        /* end QSORT_SWAP */
                    }
                    _jump_over:
                    _left_ptr = _lo + 1;
                    _right_ptr = _hi - 1;
                    /* Here's the famous ``collapse the walls'' section of quicksort. \
                         Gotta like those tight inner loops!    They are the main reason    \
                         that this algorithm runs much faster than others. */
                    do {
                        /* while QSORT_LT (_left_ptr, _mid)) */
                        while (((*_left_ptr) < (*_mid))) {
                            ++_left_ptr;
                        }
                        /* while QSORT_LT (_mid, _right_ptr)) */
                        while (((*_mid) < (*_right_ptr))) {
                            --_right_ptr;
                        }
                        if (_left_ptr < _right_ptr) {
                            /* QSORT_SWAP (_left_ptr, _right_ptr, _hold); */
                            vtmp = *_left_ptr;
                            ktmp = keys[_left_ptr - vals];
                            *_left_ptr = *_right_ptr;
                            keys[_left_ptr - vals] = keys[_right_ptr - vals];
                            *_right_ptr = vtmp;
                            keys[_right_ptr - vals] = ktmp;
                            /* end QSORT_SWAP */
                            if (_mid == _left_ptr) {
                                _mid = _right_ptr;
                            } else if (_mid == _right_ptr) {
                                _mid = _left_ptr;
                            }
                            ++_left_ptr;
                            --_right_ptr;
                        } else if (_left_ptr == _right_ptr) {
                            ++_left_ptr;
                            --_right_ptr;
                            break;
                        }
                    } while (_left_ptr <= _right_ptr);

                    /* Set up pointers for next iteration.    First determine whether     \
                        left and right partitions are below the threshold size.    If so, \
                        ignore one or both.    Otherwise, push the larger partition's     \
                        bounds on the stack and continue sorting the smaller one. */
                    if (_right_ptr - _lo <= DA_QSORT_MAX_THRESH) {
                        if (_hi - _left_ptr <= DA_QSORT_MAX_THRESH) {
                            /* Ignore both small partitions. */
                            /* QSORT_POP (_lo, _hi, _top); */
                            ((--_top, (_lo = _top->_lo), (_hi = _top->_hi)));
                        } else {
                            /* Ignore small left partition. */
                            _lo = _left_ptr;
                        }
                    } else if (_hi - _left_ptr <= DA_QSORT_MAX_THRESH) {
                        /* Ignore small right partition. */
                        _hi = _right_ptr;
                    } else if (_right_ptr - _lo > _hi - _left_ptr) {
                        /* Push larger left partition indices. */
                        /* QSORT_PUSH (_top, _lo, _right_ptr); */
                        (((_top->_lo = (_lo)), (_top->_hi = (_right_ptr)), ++_top));
                        _lo = _left_ptr;
                    } else {
                        /* Push larger right partition indices. */
                        /* QSORT_PUSH (_top, _left_ptr, _hi); */
                        (((_top->_lo = (_left_ptr)), (_top->_hi =
                                                                (_hi)), ++_top));
                        _hi = _right_ptr;
                    }
                }
            }

            /* Once the BASE array is partially sorted by quicksort the rest    \
             is completely sorted using insertion sort, since this is efficient \
             for partitions below MAX_THRESH size. BASE points to the             \
             beginning of the array to sort, and END_PTR points at the very     \
             last element in the array (*not* one beyond it!). */
            {
                V *const _end_ptr = _base + _elems - 1;
                V *_tmp_ptr = _base;
                V *_run_ptr;
                V *_thresh;
                _thresh = _base + DA_QSORT_MAX_THRESH;
                if (_thresh > _end_ptr) {
                    _thresh = _end_ptr;
                }

                /* Find smallest element in first threshold and place it at the \
                     array's beginning.    This is the smallest array element,        \
                     and the operation speeds up insertion sort's inner loop. */
                for (_run_ptr = _tmp_ptr + 1; _run_ptr <= _thresh; ++_run_ptr) {
                    if (((*_run_ptr) < (*_tmp_ptr))) { /* if QSORT_LT (_run_ptr, _tmp_ptr)) */
                        _tmp_ptr = _run_ptr;
                    }
                }
                if (_tmp_ptr != _base) {
                    /* QSORT_SWAP (_tmp_ptr, _base, _hold); */
                    vtmp = *_tmp_ptr;
                    ktmp = keys[_tmp_ptr - vals];
                    *_tmp_ptr = *_base;
                    keys[_tmp_ptr - vals] = keys[_base - vals];
                    *_base = vtmp;
                    keys[_base - vals] = ktmp;
                    /* end QSORT_SWAP */
                }

                /* Insertion sort, running from left-hand-side            \
                     up to right-hand-side.    */

                _run_ptr = _base + 1;
                while (++_run_ptr <= _end_ptr) {
                    _tmp_ptr = _run_ptr - 1;
                    /* while (QSORT_LT (_run_ptr, _tmp_ptr)) */
                    while (((*_run_ptr) < (*_tmp_ptr))) {
                        --_tmp_ptr;
                    }
                    ++_tmp_ptr;
                    if (_tmp_ptr != _run_ptr) {
                        V *_trav = _run_ptr + 1;
                        while (--_trav >= _run_ptr) {
                            V *_hi;
                            V *_lo;
                            vtmp = *_trav;
                            ktmp = keys[_trav - vals];
                            for (_hi = _lo = _trav; --_lo >= _tmp_ptr; _hi = _lo) {
                                *_hi = *_lo;
                                keys[_hi - vals] = keys[_lo - vals];
                            }
                            *_hi = vtmp;
                            keys[_hi - vals] = ktmp;
                        }
                    }
                }
            }
        }
    }


    /**
     * Sort key-value pairs in the keys and values arrays in decreasing order by value.
     * In-line qsort implementation.    Differs from traditional qsort() routine in that
     * it does not call a function to swap elements. Moreover, it sorts two arrays
     * simultaneously based on the values in the vals array and swaps the associated
     * key with each swapped value.
     *
     * @param n     Size of keys and vals arrays.
     * @param keys    Array of keys.
     * @param vals    Array of values.
     */
    template<class K, class V>
    void kvsortd(const size_t n, K *const keys, V *const vals) {
        {
            V *const _base = (vals);
            const size_t _elems = (n);
            K ktmp;
            V vtmp;
            if (_elems == 0) {
                return;
            }
            if (_elems > DA_QSORT_MAX_THRESH) {
                V *_lo = _base;
                V *_hi = _lo + _elems - 1;
                struct {
                    V *_hi;
                    V *_lo;
                } _stack[DA_QSORT_STACK_SIZE], *_top = _stack + 1;
                while ((_stack < _top)) {
                    V *_left_ptr;
                    V *_right_ptr;
                    /* Select median value from among LO, MID, and HI. Rearrange    \
                         LO and HI so the three values are sorted. This lowers the    \
                         probability of picking a pathological pivot value and        \
                         skips a comparison for both the LEFT_PTR and RIGHT_PTR in    \
                         the while loops. */
                    V *_mid = _lo + ((_hi - _lo) >> 1);
                    if (((*_mid) > (*_lo))) {    /* if (QSORT_GT (_mid, _lo)) */
                        /* QSORT_SWAP (_mid, _lo, _hold); */
                        vtmp = *_mid;
                        ktmp = keys[_mid - vals];
                        *_mid = *_lo;
                        keys[_mid - vals] = keys[_lo - vals];
                        *_lo = vtmp;
                        keys[_lo - vals] = ktmp;
                        /* end QSORT_SWAP */
                    }
                    if (((*_hi) > (*_mid))) {    /* if (QSORT_GT (_hi, _mid)) */
                        /* QSORT_SWAP (_mid, _hi, _hold); */
                        vtmp = *_mid;
                        ktmp = keys[_mid - vals];
                        *_mid = *_hi;
                        keys[_mid - vals] = keys[_hi - vals];
                        *_hi = vtmp;
                        keys[_hi - vals] = ktmp;
                        /* end QSORT_SWAP */
                    } else {
                        goto _jump_over;
                    }
                    if (((*_mid) > (*_lo))) {    /* if (QSORT_GT (_mid, _lo)) */
                        /* QSORT_SWAP (_mid, _lo, _hold); */
                        vtmp = *_mid;
                        ktmp = keys[_mid - vals];
                        *_mid = *_lo;
                        keys[_mid - vals] = keys[_lo - vals];
                        *_lo = vtmp;
                        keys[_lo - vals] = ktmp;
                        /* end QSORT_SWAP */
                    }
                    _jump_over:
                    _left_ptr = _lo + 1;
                    _right_ptr = _hi - 1;
                    /* Here's the famous ``collapse the walls'' section of quicksort. \
                         Gotta like those tight inner loops!    They are the main reason    \
                         that this algorithm runs much faster than others. */
                    do {
                        /* while QSORT_GT (_left_ptr, _mid)) */
                        while (((*_left_ptr) > (*_mid))) {
                            ++_left_ptr;
                        }
                        /* while QSORT_GT (_mid, _right_ptr)) */
                        while (((*_mid) > (*_right_ptr))) {
                            --_right_ptr;
                        }
                        if (_left_ptr < _right_ptr) {
                            /* QSORT_SWAP (_left_ptr, _right_ptr, _hold); */
                            vtmp = *_left_ptr;
                            ktmp = keys[_left_ptr - vals];
                            *_left_ptr = *_right_ptr;
                            keys[_left_ptr - vals] = keys[_right_ptr - vals];
                            *_right_ptr = vtmp;
                            keys[_right_ptr - vals] = ktmp;
                            /* end QSORT_SWAP */
                            if (_mid == _left_ptr) {
                                _mid = _right_ptr;
                            } else if (_mid == _right_ptr) {
                                _mid = _left_ptr;
                            }
                            ++_left_ptr;
                            --_right_ptr;
                        } else if (_left_ptr == _right_ptr) {
                            ++_left_ptr;
                            --_right_ptr;
                            break;
                        }
                    } while (_left_ptr <= _right_ptr);

                    /* Set up pointers for next iteration.    First determine whether     \
                        left and right partitions are below the threshold size.    If so, \
                        ignore one or both.    Otherwise, push the larger partition's     \
                        bounds on the stack and continue sorting the smaller one. */
                    if (_right_ptr - _lo <= DA_QSORT_MAX_THRESH) {
                        if (_hi - _left_ptr <= DA_QSORT_MAX_THRESH) {
                            /* Ignore both small partitions. */
                            /* QSORT_POP (_lo, _hi, _top); */
                            ((--_top, (_lo = _top->_lo), (_hi = _top->_hi)));
                        } else {
                            /* Ignore small left partition. */
                            _lo = _left_ptr;
                        }
                    } else if (_hi - _left_ptr <= DA_QSORT_MAX_THRESH) {
                        /* Ignore small right partition. */
                        _hi = _right_ptr;
                    } else if (_right_ptr - _lo > _hi - _left_ptr) {
                        /* Push larger left partition indices. */
                        /* QSORT_PUSH (_top, _lo, _right_ptr); */
                        (((_top->_lo = (_lo)), (_top->_hi = (_right_ptr)), ++_top));
                        _lo = _left_ptr;
                    } else {
                        /* Push larger right partition indices. */
                        /* QSORT_PUSH (_top, _left_ptr, _hi); */
                        (((_top->_lo = (_left_ptr)), (_top->_hi =
                                                                (_hi)), ++_top));
                        _hi = _right_ptr;
                    }
                }
            }

            /* Once the BASE array is partially sorted by quicksort the rest    \
             is completely sorted using insertion sort, since this is efficient \
             for partitions below MAX_THRESH size. BASE points to the             \
             beginning of the array to sort, and END_PTR points at the very     \
             last element in the array (*not* one beyond it!). */
            {
                V *const _end_ptr = _base + _elems - 1;
                V *_tmp_ptr = _base;
                V *_run_ptr;
                V *_thresh;
                _thresh = _base + DA_QSORT_MAX_THRESH;
                if (_thresh > _end_ptr) {
                    _thresh = _end_ptr;
                }

                /* Find smallest element in first threshold and place it at the \
                     array's beginning.    This is the smallest array element,        \
                     and the operation speeds up insertion sort's inner loop. */
                for (_run_ptr = _tmp_ptr + 1; _run_ptr <= _thresh; ++_run_ptr) {
                    if (((*_run_ptr) > (*_tmp_ptr))) { /* if QSORT_GT (_run_ptr, _tmp_ptr)) */
                        _tmp_ptr = _run_ptr;
                    }
                }
                if (_tmp_ptr != _base) {
                    /* QSORT_SWAP (_tmp_ptr, _base, _hold); */
                    vtmp = *_tmp_ptr;
                    ktmp = keys[_tmp_ptr - vals];
                    *_tmp_ptr = *_base;
                    keys[_tmp_ptr - vals] = keys[_base - vals];
                    *_base = vtmp;
                    keys[_base - vals] = ktmp;
                    /* end QSORT_SWAP */
                }

                /* Insertion sort, running from left-hand-side            \
                     up to right-hand-side.    */

                _run_ptr = _base + 1;
                while (++_run_ptr <= _end_ptr) {
                    _tmp_ptr = _run_ptr - 1;
                    /* while (QSORT_GT (_run_ptr, _tmp_ptr)) */
                    while (((*_run_ptr) > (*_tmp_ptr))) {
                        --_tmp_ptr;
                    }
                    ++_tmp_ptr;
                    if (_tmp_ptr != _run_ptr) {
                        V *_trav = _run_ptr + 1;
                        while (--_trav >= _run_ptr) {
                            V *_hi;
                            V *_lo;
                            vtmp = *_trav;
                            ktmp = keys[_trav - vals];
                            for (_hi = _lo = _trav; --_lo >= _tmp_ptr; _hi = _lo) {
                                *_hi = *_lo;
                                keys[_hi - vals] = keys[_lo - vals];
                            }
                            *_hi = vtmp;
                            keys[_hi - vals] = ktmp;
                        }
                    }
                }
            }
        }
    }




    /**
     * Sort values array in ascending order.
     * In-line qsort implementation.    Differs from traditional qsort() routine in that
     * it does not call a function to swap elements. 
     *
     * @param n     Size of keys and vals arrays.
     * @param vals    Array of values.
     */
    /* Copyright (C) 1991, 1992, 1996, 1997, 1999 Free Software Foundation, Inc.
        This file is part of the GNU C Library.
        Written by Douglas C. Schmidt (schmidt@ics.uci.edu).

        The GNU C Library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Lesser General Public
        License as published by the Free Software Foundation; either
        version 2.1 of the License, or (at your option) any later version.

        The GNU C Library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the GNU
        Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public
        License along with the GNU C Library; if not, write to the Free
        Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
        02111-1307 USA.    */
    template<class V>
    void sorti(const size_t n, V *const vals) {
        {
            V *const _base = (vals);
            const size_t _elems = (n);
            V vtmp;
            if (_elems == 0) {
                return;
            }
            if (_elems > DA_QSORT_MAX_THRESH) {
                V *_lo = _base;
                V *_hi = _lo + _elems - 1;
                struct {
                    V *_hi;
                    V *_lo;
                } _stack[DA_QSORT_STACK_SIZE], *_top = _stack + 1;
                while ((_stack < _top)) {
                    V *_left_ptr;
                    V *_right_ptr;
                    /* Select median value from among LO, MID, and HI. Rearrange    \
                         LO and HI so the three values are sorted. This lowers the    \
                         probability of picking a pathological pivot value and        \
                         skips a comparison for both the LEFT_PTR and RIGHT_PTR in    \
                         the while loops. */
                    V *_mid = _lo + ((_hi - _lo) >> 1);
                    if (((*_mid) < (*_lo))) {    /* if (QSORT_LT (_mid, _lo)) */
                        /* QSORT_SWAP (_mid, _lo, _hold); */
                        vtmp = *_mid;
                        *_mid = *_lo;
                        *_lo = vtmp;
                        /* end QSORT_SWAP */
                    }
                    if (((*_hi) < (*_mid))) {    /* if (QSORT_LT (_hi, _mid)) */
                        /* QSORT_SWAP (_mid, _hi, _hold); */
                        vtmp = *_mid;
                        *_mid = *_hi;
                        *_hi = vtmp;
                        /* end QSORT_SWAP */
                    } else {
                        goto _jump_over;
                    }
                    if (((*_mid) < (*_lo))) {    /* if (QSORT_LT (_mid, _lo)) */
                        /* QSORT_SWAP (_mid, _lo, _hold); */
                        vtmp = *_mid;
                        *_mid = *_lo;
                        *_lo = vtmp;
                        /* end QSORT_SWAP */
                    }
                    _jump_over:
                    _left_ptr = _lo + 1;
                    _right_ptr = _hi - 1;
                    /* Here's the famous ``collapse the walls'' section of quicksort. \
                         Gotta like those tight inner loops!    They are the main reason    \
                         that this algorithm runs much faster than others. */
                    do {
                        /* while QSORT_LT (_left_ptr, _mid)) */
                        while (((*_left_ptr) < (*_mid))) {
                            ++_left_ptr;
                        }
                        /* while QSORT_LT (_mid, _right_ptr)) */
                        while (((*_mid) < (*_right_ptr))) {
                            --_right_ptr;
                        }
                        if (_left_ptr < _right_ptr) {
                            /* QSORT_SWAP (_left_ptr, _right_ptr, _hold); */
                            vtmp = *_left_ptr;
                            *_left_ptr = *_right_ptr;
                            *_right_ptr = vtmp;
                            /* end QSORT_SWAP */
                            if (_mid == _left_ptr) {
                                _mid = _right_ptr;
                            } else if (_mid == _right_ptr) {
                                _mid = _left_ptr;
                            }
                            ++_left_ptr;
                            --_right_ptr;
                        } else if (_left_ptr == _right_ptr) {
                            ++_left_ptr;
                            --_right_ptr;
                            break;
                        }
                    } while (_left_ptr <= _right_ptr);

                    /* Set up pointers for next iteration.    First determine whether     \
                        left and right partitions are below the threshold size.    If so, \
                        ignore one or both.    Otherwise, push the larger partition's     \
                        bounds on the stack and continue sorting the smaller one. */
                    if (_right_ptr - _lo <= DA_QSORT_MAX_THRESH) {
                        if (_hi - _left_ptr <= DA_QSORT_MAX_THRESH) {
                            /* Ignore both small partitions. */
                            /* QSORT_POP (_lo, _hi, _top); */
                            ((--_top, (_lo = _top->_lo), (_hi = _top->_hi)));
                        } else {
                            /* Ignore small left partition. */
                            _lo = _left_ptr;
                        }
                    } else if (_hi - _left_ptr <= DA_QSORT_MAX_THRESH) {
                        /* Ignore small right partition. */
                        _hi = _right_ptr;
                    } else if (_right_ptr - _lo > _hi - _left_ptr) {
                        /* Push larger left partition indices. */
                        /* QSORT_PUSH (_top, _lo, _right_ptr); */
                        (((_top->_lo = (_lo)), (_top->_hi = (_right_ptr)), ++_top));
                        _lo = _left_ptr;
                    } else {
                        /* Push larger right partition indices. */
                        /* QSORT_PUSH (_top, _left_ptr, _hi); */
                        (((_top->_lo = (_left_ptr)), (_top->_hi = (_hi)), ++_top));
                        _hi = _right_ptr;
                    }
                }
            }

            /* Once the BASE array is partially sorted by quicksort the rest    \
             is completely sorted using insertion sort, since this is efficient \
             for partitions below MAX_THRESH size. BASE points to the             \
             beginning of the array to sort, and END_PTR points at the very     \
             last element in the array (*not* one beyond it!). */
            {
                V *const _end_ptr = _base + _elems - 1;
                V *_tmp_ptr = _base;
                V *_run_ptr;
                V *_thresh;
                _thresh = _base + DA_QSORT_MAX_THRESH;
                if (_thresh > _end_ptr) {
                    _thresh = _end_ptr;
                }

                /* Find smallest element in first threshold and place it at the \
                     array's beginning.    This is the smallest array element,        \
                     and the operation speeds up insertion sort's inner loop. */
                for (_run_ptr = _tmp_ptr + 1; _run_ptr <= _thresh; ++_run_ptr) {
                    if (((*_run_ptr) < (*_tmp_ptr))) { /* if QSORT_LT (_run_ptr, _tmp_ptr)) */
                        _tmp_ptr = _run_ptr;
                    }
                }
                if (_tmp_ptr != _base) {
                    /* QSORT_SWAP (_tmp_ptr, _base, _hold); */
                    vtmp = *_tmp_ptr;
                    *_tmp_ptr = *_base;
                    *_base = vtmp;
                    /* end QSORT_SWAP */
                }

                /* Insertion sort, running from left-hand-side            \
                     up to right-hand-side.    */

                _run_ptr = _base + 1;
                while (++_run_ptr <= _end_ptr) {
                    _tmp_ptr = _run_ptr - 1;
                    /* while (QSORT_LT (_run_ptr, _tmp_ptr)) */
                    while (((*_run_ptr) < (*_tmp_ptr))) {
                        --_tmp_ptr;
                    }
                    ++_tmp_ptr;
                    if (_tmp_ptr != _run_ptr) {
                        V *_trav = _run_ptr + 1;
                        while (--_trav >= _run_ptr) {
                            V *_hi;
                            V *_lo;
                            vtmp = *_trav;
                            for (_hi = _lo = _trav; --_lo >= _tmp_ptr; _hi = _lo) {
                                *_hi = *_lo;
                            }
                            *_hi = vtmp;
                        }
                    }
                }
            }
        }
    }


    /**
     * Sort values array in decreasing order by value.
     * In-line qsort implementation.    Differs from traditional qsort() routine in that
     * it does not call a function to swap elements. 
     *
     * @param n     Size of keys and vals arrays.
     * @param keys    Array of keys.
     * @param vals    Array of values.
     */
    template<class V>
    void sortd(const size_t n, V *const vals) {
        {
            V *const _base = (vals);
            const size_t _elems = (n);
            V vtmp;
            if (_elems == 0) {
                return;
            }
            if (_elems > DA_QSORT_MAX_THRESH) {
                V *_lo = _base;
                V *_hi = _lo + _elems - 1;
                struct {
                    V *_hi;
                    V *_lo;
                } _stack[DA_QSORT_STACK_SIZE], *_top = _stack + 1;
                while ((_stack < _top)) {
                    V *_left_ptr;
                    V *_right_ptr;
                    /* Select median value from among LO, MID, and HI. Rearrange    \
                         LO and HI so the three values are sorted. This lowers the    \
                         probability of picking a pathological pivot value and        \
                         skips a comparison for both the LEFT_PTR and RIGHT_PTR in    \
                         the while loops. */
                    V *_mid = _lo + ((_hi - _lo) >> 1);
                    if (((*_mid) > (*_lo))) {    /* if (QSORT_GT (_mid, _lo)) */
                        /* QSORT_SWAP (_mid, _lo, _hold); */
                        vtmp = *_mid;
                        *_mid = *_lo;
                        *_lo = vtmp;
                        /* end QSORT_SWAP */
                    }
                    if (((*_hi) > (*_mid))) {    /* if (QSORT_GT (_hi, _mid)) */
                        /* QSORT_SWAP (_mid, _hi, _hold); */
                        vtmp = *_mid;
                        *_mid = *_hi;
                        *_hi = vtmp;
                        /* end QSORT_SWAP */
                    } else {
                        goto _jump_over;
                    }
                    if (((*_mid) > (*_lo))) {    /* if (QSORT_GT (_mid, _lo)) */
                        /* QSORT_SWAP (_mid, _lo, _hold); */
                        vtmp = *_mid;
                        *_mid = *_lo;
                        *_lo = vtmp;
                        /* end QSORT_SWAP */
                    }
                    _jump_over:
                    _left_ptr = _lo + 1;
                    _right_ptr = _hi - 1;
                    /* Here's the famous ``collapse the walls'' section of quicksort. \
                         Gotta like those tight inner loops!    They are the main reason    \
                         that this algorithm runs much faster than others. */
                    do {
                        /* while QSORT_GT (_left_ptr, _mid)) */
                        while (((*_left_ptr) > (*_mid))) {
                            ++_left_ptr;
                        }
                        /* while QSORT_GT (_mid, _right_ptr)) */
                        while (((*_mid) > (*_right_ptr))) {
                            --_right_ptr;
                        }
                        if (_left_ptr < _right_ptr) {
                            /* QSORT_SWAP (_left_ptr, _right_ptr, _hold); */
                            vtmp = *_left_ptr;
                            *_left_ptr = *_right_ptr;
                            *_right_ptr = vtmp;
                            /* end QSORT_SWAP */
                            if (_mid == _left_ptr) {
                                _mid = _right_ptr;
                            } else if (_mid == _right_ptr) {
                                _mid = _left_ptr;
                            }
                            ++_left_ptr;
                            --_right_ptr;
                        } else if (_left_ptr == _right_ptr) {
                            ++_left_ptr;
                            --_right_ptr;
                            break;
                        }
                    } while (_left_ptr <= _right_ptr);

                    /* Set up pointers for next iteration.    First determine whether     \
                        left and right partitions are below the threshold size.    If so, \
                        ignore one or both.    Otherwise, push the larger partition's     \
                        bounds on the stack and continue sorting the smaller one. */
                    if (_right_ptr - _lo <= DA_QSORT_MAX_THRESH) {
                        if (_hi - _left_ptr <= DA_QSORT_MAX_THRESH) {
                            /* Ignore both small partitions. */
                            /* QSORT_POP (_lo, _hi, _top); */
                            ((--_top, (_lo = _top->_lo), (_hi = _top->_hi)));
                        } else {
                            /* Ignore small left partition. */
                            _lo = _left_ptr;
                        }
                    } else if (_hi - _left_ptr <= DA_QSORT_MAX_THRESH) {
                        /* Ignore small right partition. */
                        _hi = _right_ptr;
                    } else if (_right_ptr - _lo > _hi - _left_ptr) {
                        /* Push larger left partition indices. */
                        /* QSORT_PUSH (_top, _lo, _right_ptr); */
                        (((_top->_lo = (_lo)), (_top->_hi = (_right_ptr)), ++_top));
                        _lo = _left_ptr;
                    } else {
                        /* Push larger right partition indices. */
                        /* QSORT_PUSH (_top, _left_ptr, _hi); */
                        (((_top->_lo = (_left_ptr)), (_top->_hi = (_hi)), ++_top));
                        _hi = _right_ptr;
                    }
                }
            }

            /* Once the BASE array is partially sorted by quicksort the rest    \
             is completely sorted using insertion sort, since this is efficient \
             for partitions below MAX_THRESH size. BASE points to the             \
             beginning of the array to sort, and END_PTR points at the very     \
             last element in the array (*not* one beyond it!). */
            {
                V *const _end_ptr = _base + _elems - 1;
                V *_tmp_ptr = _base;
                V *_run_ptr;
                V *_thresh;
                _thresh = _base + DA_QSORT_MAX_THRESH;
                if (_thresh > _end_ptr) {
                    _thresh = _end_ptr;
                }

                /* Find smallest element in first threshold and place it at the \
                     array's beginning.    This is the smallest array element,        \
                     and the operation speeds up insertion sort's inner loop. */
                for (_run_ptr = _tmp_ptr + 1; _run_ptr <= _thresh; ++_run_ptr) {
                    if (((*_run_ptr) > (*_tmp_ptr))) { /* if QSORT_GT (_run_ptr, _tmp_ptr)) */
                        _tmp_ptr = _run_ptr;
                    }
                }
                if (_tmp_ptr != _base) {
                    /* QSORT_SWAP (_tmp_ptr, _base, _hold); */
                    vtmp = *_tmp_ptr;
                    *_tmp_ptr = *_base;
                    *_base = vtmp;
                    /* end QSORT_SWAP */
                }

                /* Insertion sort, running from left-hand-side            \
                     up to right-hand-side.    */

                _run_ptr = _base + 1;
                while (++_run_ptr <= _end_ptr) {
                    _tmp_ptr = _run_ptr - 1;
                    /* while (QSORT_GT (_run_ptr, _tmp_ptr)) */
                    while (((*_run_ptr) > (*_tmp_ptr))) {
                        --_tmp_ptr;
                    }
                    ++_tmp_ptr;
                    if (_tmp_ptr != _run_ptr) {
                        V *_trav = _run_ptr + 1;
                        while (--_trav >= _run_ptr) {
                            V *_hi;
                            V *_lo;
                            vtmp = *_trav;
                            for (_hi = _lo = _trav; --_lo >= _tmp_ptr; _hi = _lo) {
                                *_hi = *_lo;
                            }
                            *_hi = vtmp;
                        }
                    }
                }
            }
        }
    }

    /**
     * Select-K routines.
     */

    #define VSWAP(a, b, stmp) do { stmp = (a); (a) = (b); (b) = stmp; } while (0)

    /**
     * This function puts the `k` largest values in the beginning of the array
     * @param n      Size of array to select from
     * @param varr   Array of values to select from
     * @param k      Number of smallest values to move to the start of the vector.
     */
    template <typename V>
    void selectd(
        const size_t n,
        V* varr,
        const size_t k
    ){
        size_t i, j, lo, hi, mid;
        V vtmp;
        V pivot;

        if (n <= k){
            return;         /* return if the vector has fewer elements than we want */
        }

        for (lo = 0, hi = n - 1; lo < hi;)
        {
            mid = lo + ((hi - lo) >> 1);

            /* select the median */
            if (varr[lo] < varr[mid]){
                mid = lo;
            }
            if (varr[hi] > varr[mid]){
                mid = hi;
            } else{
                goto jump_over;
            }
            if (varr[lo] < varr[mid]){
                mid = lo;
            }

        jump_over:
            //swap mid and high
            VSWAP(varr[mid], varr[hi], vtmp);
            pivot = varr[hi];

            /* the partitioning algorithm */
            for (i = lo - 1, j = lo; j < hi; j++)
            {
                if (varr[j] >= pivot)
                {
                    i++;
                    VSWAP(varr[i], varr[j], vtmp);
                }
            }
            i++;
            VSWAP(varr[i], varr[hi], vtmp);

            if (i > k){
                hi = i - 1;
            } else if (i < k){
                lo = i + 1;
            } else{
                break;
            }
        }
    }

    /**
     * This function puts the `k` smallest values in the beginning of the vector.
     * @param n      Size of array to select from
     * @param varr   Array of values to select from
     * @param k      Number of smallest values to move to the start of the vector.
     */
    template <typename V>
    void selecti(
        const size_t n,
        V* varr,
        const size_t k
    ){
        size_t i, j, lo, hi, mid;
        V vtmp;
        V pivot;

        if (n <= k){
            return;         /* return if the vector has fewer elements than we want */
        }

        for (lo = 0, hi = n - 1; lo < hi;)
        {
            mid = lo + ((hi - lo) >> 1);

            /* select the median */
            if (varr[lo] > varr[mid]){
                mid = lo;
            }
            if (varr[hi] < varr[mid]){
                mid = hi;
            } else{
                goto jump_over;
            }
            if (varr[lo] > varr[mid]){
                mid = lo;
            }

        jump_over:
            //swap mid and high
            VSWAP(varr[mid], varr[hi], vtmp);
            pivot = varr[hi];

            /* the partitioning algorithm */
            for (i = lo - 1, j = lo; j < hi; j++)
            {
                if (varr[j] <= pivot)
                {
                    i++;
                    VSWAP(varr[i], varr[j], vtmp);
                }
            }
            i++;
            VSWAP(varr[i], varr[hi], vtmp);

            if (i > k){
                hi = i - 1;
            } else if (i < k){
                lo = i + 1;
            } else{
                break;
            }
        }
    }

    #undef VSWAP
  
    /**
     * Sort values array in ascending/decreasing order by value.
     * In-line qsort implementation.    Differs from traditional qsort() routine in that
     * it does not call a function to swap elements. 
     *
     * @param n     Size of keys and vals arrays.
     * @param keys    Array of keys.
     * @param vals    Array of values.
     * @param ascending Whether to select the smallest or highest values
     */
    template<class V>
    void sort(const size_t n, V *const vals, const bool ascending) {
        if(ascending){
            sorti(n, vals);
        } else {
            sortd(n, vals);
        }
    }


    /**
     * This function puts the `k` smallest/highest values in the beginning of the vector.
     * @param n      Size of array to select from
     * @param varr   Array of values to select from
     * @param k      Number of smallest values to move to the start of the vector.
     * @param ascending Whether to select the smallest or highest values
     */
    template <typename V>
    void select(
        const size_t n,
        V* varr,
        const size_t k,
        bool ascending
    ){
        if(ascending){
            selecti(n, varr, k);
        } else {
            selectd(n, varr, k);
        }
    }


    /**
     * This function puts the `k` largest values in the beginning of the array and sorts them.
     * @param n      Size of array to select from
     * @param varr   Array of values to select from
     * @param k      Number of smallest values to move to the start of the vector.
     */
    template <typename V>
    void selectsortd(
        const size_t n,
        V* varr,
        const size_t k
    ){
        selectd(n, varr, k);
        sortd(std::min(k, n), varr);
    }
  
    /**
     * This function puts the `k` smallest values in the beginning of the array and sorts them.
     * @param n      Size of array to select from
     * @param varr   Array of values to select from
     * @param k      Number of smallest values to move to the start of the vector.
     */
    template <typename V>
    void selectsorti(
        const size_t n,
        V* varr,
        const size_t k
    ){
        selecti(n, varr, k);
        sorti(std::min(k, n), varr);
    }
  
    /**
     * This function puts the `k` smallest/highest values in the beginning of the array and sorts them.
     * @param n      Size of array to select from
     * @param varr   Array of values to select from
     * @param k      Number of smallest values to move to the start of the vector.
     * @param ascending Whether to select the smallest or highest values
     */
    template <typename V>
    void selectsort(
        const size_t n,
        V* varr,
        const size_t k,
        bool ascending
    ){
        if(ascending){
            selecti(n, varr, k);
            sorti(std::min(k, n), varr);
        } else {
            selectd(n, varr, k);
            sortd(std::min(k, n), varr);
        }
    }


} /* end namespace */
