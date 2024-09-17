#ifndef _NG_MIXED_BUFFER_HPP_
#define _NG_MIXED_BUFFER_HPP_
#include <vector>
#include "ng_typedefs.hpp"
namespace ngineer
{
    constexpr char MIXED_BUFFER_STACK_SIZE = 16; // enough for a 4x4 matrix

    /// @brief  A buffer that stores some memory on the heap before 
    ///         resorting to using memory on the stack. This amount
    ///         is determined by `MIXED_BUFFER_STACK_SIZE`, and (if
    ///         edited) should reflect a number of elements that 
    ///         covers common cases 
    /// @tparam T The type
    template<class T>
    class MixedBuffer
    {
    private:
        buflen_t length;
        buflen_t capacity;
        T stack_buffer[MIXED_BUFFER_STACK_SIZE];
        std::vector<T> heap_buffer;
    public:
        MixedBuffer();
        MixedBuffer(std::vector<T> elements);
        MixedBuffer(buflen_t capacity);
        ~MixedBuffer();
        inline buflen_t get_length(void) const noexcept;
        void push_back(T element);
        T& back(void);
        void pop_back(void);
        T& operator[](buflen_t index);

        // NgIterator implementations
        T *begin(void);
        T *end(void);
    };

    template<class T>
    MixedBuffer<T>::MixedBuffer()
    {
        length = 0;
        capacity = MIXED_BUFFER_STACK_SIZE;
        heap_buffer = std::vector<T>();
        for (buflen_t i = 0; i < MIXED_BUFFER_STACK_SIZE; i++)
        {
            stack_buffer[i] = 0x00;
        }
    }

    template<class T>
    MixedBuffer<T>::MixedBuffer(std::vector<T> elements)
    {   
        MixedBuffer(); // Initialize to defaults
        while(0 != elements.size())
        {
            push_back(elements.back());
            elements.pop_back();
        }
    }


    template<class T>
    MixedBuffer<T>::MixedBuffer(buflen_t cap)
    {
        capacity = cap;
        length = 0;

        for (buflen_t i = 0; i < MIXED_BUFFER_STACK_SIZE; i++)
        {
            stack_buffer[i] = 0x00;
        }
        if (capacity >= MIXED_BUFFER_STACK_SIZE)
        {
            heap_buffer.reserve(capacity - MIXED_BUFFER_STACK_SIZE);
        }
    }

    template<class T>
    MixedBuffer<T>::~MixedBuffer()
    {
    }

    template<class T>
    inline buflen_t MixedBuffer<T>::get_length(void) const noexcept
    {
        return length;
    }

    template<class T>
    void MixedBuffer<T>::push_back(T element)
    {
        if (length >= MIXED_BUFFER_STACK_SIZE)
        {
            heap_buffer.push_back(element);
        }
        else
        {
            stack_buffer[length] = element;
        }
        length++;
    }

    template<class T>
    T& MixedBuffer<T>::back(void)
    {
        return this[length - 1];
    }

    template<class T>
    void MixedBuffer<T>::pop_back(void)
    {
        if (length >= MIXED_BUFFER_STACK_SIZE)
        {
            heap_buffer.pop_back();
        }
        length--;
    }

    template<class T>
    T& MixedBuffer<T>::operator[](buflen_t index)
    {
        if (length >= MIXED_BUFFER_STACK_SIZE)
        {
            return heap_buffer[index - MIXED_BUFFER_STACK_SIZE];
        }
        else
        {
            return stack_buffer[index];
        }
    }

    template<class T>
    T *MixedBuffer<T>::begin() 
    {
        return stack_buffer;
    }

    template<class T>
    T *MixedBuffer<T>::end()
    {
        return heap_buffer.end();
    }
} // namespace ngineer
#endif // _MIXED_BUFFER_HPP_