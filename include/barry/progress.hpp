#ifndef BARRY_PROGRESS_HPP
#define BARRY_PROGRESS_HPP

#ifndef BARRY_PROGRESS_BAR_WIDTH
#define BARRY_PROGRESS_BAR_WIDTH 80
#endif

/**
 * @brief A simple progress bar
  */
class Progress {
private:
    int    width;     ///< Total width size (number of bars)
    int    n;         ///< Total number of iterations
    double step_size; ///< Size of the step
    int last_loc;     ///< Last location of the bar
    int cur_loc;      ///< Last location of the bar
    int i;            ///< Current iteration step
    
public:

    Progress(int n_, int width_);
    ~Progress() {};

    void next();
    void end();

};

inline Progress::Progress(int n_, int width_) {


    width     = std::max(7, width_ - 7);
    n         = n_;
    step_size = static_cast<double>(width)/static_cast<double>(n);
    last_loc  = 0;
    i         = 0;

}

inline void Progress::next() {

    cur_loc = std::floor((++i) * step_size);

    for (int j = 0; j < (cur_loc - last_loc); ++j)
        printf_barry("|");

    last_loc = cur_loc;

}

inline void Progress::end() {

    printf_barry(" done.\n");

}

#endif