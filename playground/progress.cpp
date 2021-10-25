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
            printf_barry("|");
    }

};

inline void Progress::end() {

    int reminder = static_cast<int>(width) - n_bars * n_steps;
    for (int j = 0; j < reminder; ++j)
        printf_barry("|");
    
    printf_barry(" done.\n");

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

    std::cout << "Length 121\n";
    Progress p5(100, 121);
    for (int i = 0; i < 100; ++i) 
        p5.next();
    
    p5.end();

    std::cout << "Checking different sizes" << std::endl;

    std::cout << "Length 20\n";
    Progress p_(20, 80);
    for (int i = 0; i < 20; ++i) 
        p_.next();
    
    p_.end();

    std::cout << "Length 50\n";
    Progress p2_(50, 80);
    for (int i = 0; i < 50; ++i) 
        p2_.next();
    
    p2_.end();

    std::cout << "Length 100\n";
    Progress p3_(100, 80);
    for (int i = 0; i < 100; ++i) 
        p3_.next();
    
    p3_.end();

    std::cout << "Length 80\n";
    Progress p4_(80, 80);
    for (int i = 0; i < 80; ++i) 
        p4_.next();
    
    p4_.end();

    std::cout << "Length 121\n";
    Progress p5_(121, 80);
    for (int i = 0; i < 121; ++i) 
        p5_.next();
    
    p5_.end();

    return 0;

}
