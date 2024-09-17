#ifndef _NG_TYPEDEFS_HPP_
#define _NG_TYPEDEFS_HPP_

namespace ngineer
{
    /// @brief A type for storing the lenghths of various buffers.
    ///        This defaults to `unsigned int` for a balance of
    ///        allowing reasonably large buffers while not making
    ///        the size of the buffer object itself larger than it
    ///        needs to be.
    typedef unsigned int buflen_t;

    /// @brief An interface for 
    /// @tparam T 
    template<class T>
    class NgIterator
    {
    protected:
        T *iterator_item;
        void iterator_increment_routine(void) {};
    public:
        NgIterator(): iterator_item(nullptr) {};
        T *operator++();
    }; // class NgIteratorInterface

    template<class T>
    T *NgIterator<T>::operator++()
    {
        iterator_increment_routine();
        return iterator_item;
    }

} // namespace ngineer

#endif // _NG_TYPEDEFS_HPP_