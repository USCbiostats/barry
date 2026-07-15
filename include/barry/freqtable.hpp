#ifndef BARRY_STATSDB_HPP 
#define BARRY_STATSDB_HPP 1
  
/**
 * @brief Frequency table of vectors
 * 
 * This is mostly used in `Support`. The main data is contained in the
 * `data` double vector. The matrix is stored in a row-wise fashion, where
 * the first element is the frequency with which the vector is observed.
 * 
 * For example, in a model with `k` terms the first k + 1 elements of
 * `data` would be:
 * 
 * - weights
 * - term 1
 * - term 2
 * - ...
 * - term k
 * 
 */
template<typename T = double>
class FreqTable {
private:

    std::unordered_map<size_t, size_t> index;
    std::vector< double > data;
    size_t k = 0u;
    size_t n = 0u;

    /**
     * Rows whose hash collided with a row already in `index`
     * (same 64-bit hash, different statistics). Keyed by the shared
     * hash, holding the positions of the extra rows in `data`.
     * Empty unless a true hash collision occurs, so it adds no
     * per-row overhead in the common case.
     */
    std::unordered_map<size_t, std::vector<size_t>> collided;

    typename std::unordered_map<size_t, size_t>::iterator iter;

    /**
     * Whether the row stored at `pos` (position of its weight in
     * `data`) holds exactly the statistics in `x`.
     */
    bool row_matches(size_t pos, const std::vector< T > & x) const
    {
        for (size_t j = 0u; j < k; ++j)
            if (data[pos + 1u + j] != x[j])
                return false;
        return true;
    }

    /**
     * Fallback for a true hash collision: the row at `pos` shares the
     * hash `h` with `x` but holds different statistics. Looks `x` up in
     * (or appends it to) the per-hash overflow chain.
     */
    size_t add_collided(size_t h, const std::vector< T > & x)
    {

        auto & chain = collided[h];
        for (auto p : chain)
        {
            if (row_matches(p, x))
            {
                data[p] += 1.0;
                return h;
            }
        }

        chain.push_back(data.size());
        data.push_back(1.0);
        data.insert(data.end(), x.begin(), x.end());
        n++;

        return h;

    }

public:
    // size_t ncols;
    FreqTable() {};
    ~FreqTable() {};
    
    size_t add(const std::vector< T > & x, size_t * h_precomp);
    
    Counts_type                 as_vector() const;
    const std::vector< double > & get_data() const {return data;};
    const std::unordered_map<size_t,size_t> & get_index() const {return index;};
    
    void clear();
    void reserve(size_t n, size_t k);
    void print() const;

    /**
     * @brief Number of unique elements in the table.
     * (
     * @return size_t 
     */
    size_t size() const noexcept;

    size_t make_hash(const std::vector< T > & x) const;
    
};

template<typename T>  
inline size_t FreqTable<T>::add(
    const std::vector< T > & x,
    size_t * h_precomp
    ) { 
    
    // The term exists, then we add it to the list and we initialize it
    // with a single count
    size_t h;
    if (h_precomp == nullptr)
        h = make_hash(x);
    else
        h = *h_precomp;

    if (k == 0u)
    {

        index.insert({h, 0u});

        data.push_back(1.0);
        data.insert(data.end(), x.begin(), x.end());

        k = x.size();
        n++;

        return h;

    }
    else
    {


        if (x.size() != k)
            throw std::length_error(
                "The value you are trying to add doesn't have the same length used in the database."
                );

        #if __cplusplus > 201700L
        auto iter2 = index.try_emplace(h, data.size());

        if (!iter2.second)
        {

            // The hash existed: make sure the row it points to actually
            // holds the same statistics before counting it. Otherwise it
            // is a hash collision and the row goes to the overflow chain.
            size_t pos = (iter2.first)->second;
            if (row_matches(pos, x))
                data[pos] += 1.0;
            else
                return add_collided(h, x);

        }
        else
        {
            data.push_back(1.0);
            data.insert(data.end(), x.begin(), x.end());
            n++;
        }
        #else
        iter = index.find(h);

        if (iter == index.end())
        {


            index.insert({h, data.size()});
            data.push_back(1.0);
            data.insert(data.end(), x.begin(), x.end());

            n++;

            return h;

        }

        // The hash existed: make sure the row it points to actually
        // holds the same statistics before counting it. Otherwise it
        // is a hash collision and the row goes to the overflow chain.
        if (row_matches((*iter).second, x))
            data[(*iter).second] += 1.0;
        else
            return add_collided(h, x);

        #endif


    }

    return h;

}

template<typename T>
inline Counts_type FreqTable<T>::as_vector() const
{ 
    
    Counts_type ans;

    ans.reserve(index.size());

    for (size_t i = 0u; i < n; ++i)
    {
        
        std::vector< double > tmp(k, 0.0);

        for (size_t j = 1u; j < (k + 1u); ++j)
            tmp[j - 1u] = data[i * (k + 1) + j];
        
        ans.push_back(
            std::make_pair<std::vector<double>,size_t>(
                std::move(tmp),
                static_cast<size_t>(data[i * (k + 1u)])
                )
        );

    }
    
    
    return ans;
}

template<typename T>
inline void FreqTable<T>::clear()
{

    index.clear();
    collided.clear();
    data.clear();

    n = 0u;
    k = 0u;

    return;

}

template<typename T>
inline void FreqTable<T>::reserve(
    size_t n,
    size_t k
)
{

    // Figuring out the max size
    auto nk = std::min(BARRY_MAX_NUM_ELEMENTS, n * k);
    n = nk / k;
    data.reserve(nk);
    index.reserve(n);

    return;

}

// inline void StatsDB::rehash() {
//   stats.rehash();
//   return;
// }

template<typename T>
inline void FreqTable<T>::print() const
{

    size_t grand_total = 0u;

    printf_barry("%7s | %s\n", "Counts", "Stats");

    for (size_t i = 0u; i < n; ++i)
    {

        printf_barry("%7i | ", static_cast<int>(data[i * (k + 1u)]));

        for (size_t j = 1u; j < (k + 1u); ++j)
        {
            printf_barry(" %.2f", data[i * (k + 1) + j]);
        }
        printf_barry("\n");

        grand_total += static_cast<size_t>(data[i * (k + 1u)]);

    }

    printf_barry("Grand total: %i\n", static_cast<int>(grand_total));

    return;

}

template<typename T>
inline size_t FreqTable<T>::size() const noexcept
{

    // Not index.size(): rows stored in the collision overflow chain
    // have no entry of their own in -index-.
    return n;

}

template<typename T>
inline size_t FreqTable<T>::make_hash(const std::vector< T > & x) const
{

    // Well-defined for empty input (the loop below reads x[0u]).
    if (x.empty())
        return 0u;

    std::hash< T > hasher;

    // Finalizer from splitmix64: std::hash of arithmetic types is usually
    // the identity (or a bit-cast for doubles, where small integers leave
    // the low ~50 bits zero). Without this avalanche step, the low bits of
    // the combined hash carry almost no entropy, and unordered containers
    // that mask the hash with a power-of-two bucket count put (nearly) all
    // keys in a single bucket, turning lookups into O(n). The mixing runs
    // in an explicit 64-bit accumulator so the avalanche is identical where
    // size_t is 32-bit (the constants and shifts are 64-bit).
    auto mix64 = [](std::uint64_t z) -> std::uint64_t {
        z += 0x9e3779b97f4a7c15ULL;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    };

    std::uint64_t hash = mix64(hasher(x[0u]));

    // ^ makes bitwise XOR
    // 0x9e3779b9 is a 32 bit constant (comes from the golden ratio)
    // << is a shift operator, something like lhs * 2^(rhs)
    if (x.size() > 1u)
        for (size_t i = 1u; i < x.size(); ++i)
            hash ^= mix64(hasher(x[i])) + 0x9e3779b9 + (hash<<6) + (hash>>2);

    return static_cast< size_t >(hash);

}

#endif
