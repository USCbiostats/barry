#include<iostream>
#include<string>
#include<cmath>

class Progress {
private:
    int width;
    int n;
    int step;
    int n_steps;
    int step_size;
    int n_bars;
public:

    Progress(int n_, int width_);
    ~Progress() {};

    void next();
    void end();

};

inline Progress::Progress(int n_, int width_) {


    width   = width_ - 7;
    n       = n_;
    step    = 0;
    n_steps = (n > width) ?
        static_cast<int>(floor(static_cast<double>(n) / width)) : n;

    step_size  = static_cast<int>(n / n_steps);

    n_bars = static_cast<int>(std::max(1.0, floor(width / static_cast<double>(n_steps))));
};

inline void Progress::next() {

    if (!(step++ % step_size))
    {
        for (int j = 0; j < n_bars; ++j)
            printf("|");
    }

};

inline void Progress::end() {

    int reminder = static_cast<int>(width) - n_bars * n_steps;
    for (int j = 0; j < reminder; ++j)
        printf("|");
    
    printf(" done.\n");

};

int main() {

    std::cout << "Length 20\n";
    Progress p(100, 20);
    for (int i = 0; i < 100; ++i) 
        p.next();
    
    p.end();

    std::cout << "Length 50\n";
    Progress p2(100, 50);
    for (int i = 0; i < 100; ++i) 
        p2.next();
    
    p2.end();

    std::cout << "Length 100\n";
    Progress p3(100, 100);
    for (int i = 0; i < 100; ++i) 
        p3.next();
    
    p3.end();

    std::cout << "Length 80\n";
    Progress p4(100, 80);
    for (int i = 0; i < 100; ++i) 
        p4.next();
    
    p4.end();

    return 0;

}
