#ifndef __SMARTPTR_HPP__
#define __SMARTPTR_HPP__

#include "sgerrcheck.hpp"

namespace snowgoose {

    template < class T > struct SmartArrayPtr {
    public:

        /**
         * The constructor
         */
        SmartArrayPtr(unsigned int num = 0) {
            setDim(0);
            alloc(num);
        }

        /**
         * Destructor
         *
         */
        ~SmartArrayPtr() {
            dismiss();
        }

        /**
         * Copy constructor
         *
         * @param sp pointer to copy
         */
        SmartArrayPtr(const SmartArrayPtr & sp) {
            mDim = 0;
            assign(sp);
        }

        /**
         * Assignment operator
         * @param sp
         * @return
         */
        SmartArrayPtr & operator=(const SmartArrayPtr & sp) {
            dismiss();
            assign(sp);
            return *this;
        }

        /**
         * Copies the content of the array pointer
         * @param sp array to clone
         */
        void copy(const SmartArrayPtr & sp) {
            dismiss();
            alloc(sp.mDim);
            for (int i = 0; i < mDim; i++)
                mX[i] = sp[i];
        }

        /**
         * Allocate memory for this pointer
         *
         * @param num number of items
         *
         */
        void alloc(unsigned int num) {
            if (mDim != 0) {
                SG_ERROR_REPORT("Trying to allocate non-null pointer");
            }
            if (num > 0) {
                setDim(num);
                mCounter = (int*) malloc(sizeof (int));
                *((int*) mCounter) = 1;
                mX = (T*) malloc(sizeof (T) * num);
            }
        }

        /**
         * Allocate memory if not already done
         * @param num number of items
         */
        void lazyAlloc(int num) {
            if (mDim == 0)
                alloc(num);
        }

        /**
         * Reallocate memory if the storage is not sufficient
         * @param num number of items 
         */
        void realloc(int num) {
            if (mDim < num) {
                if (mDim > 0)
                    dismiss();
                alloc(num);
            }
        }

        /**
         * Pointer value access
         */
        operator T* () const {
            return (T*) mX;
        }

        /**
         *Check nullity
         */
        bool isNil() const {
            return (mDim == 0);
        }

        /**
         * Get size
         * @return size of the array
         */
        unsigned int size() const {
            return mDim;
        }
        
    private:

        void assign(const SmartArrayPtr& sp) {
            if (sp.mDim > 0) {
                mCounter = sp.mCounter;
                mX = sp.mX;
                (*(int*) (mCounter))++;
            }
            setDim(sp.mDim);
        }

        void dismiss() {
            if (mDim > 0) {
                (*(int*)mCounter)--;
                if ((*(int*) mCounter) == 0) {
                    free(mX);
                    free(mCounter);
                }
                setDim(0);
            }
        }

        void setDim(unsigned int dim) {
            mDim = dim;
        }

        int* mCounter;
        T* mX;
        unsigned int mDim;
    };

}

#endif
