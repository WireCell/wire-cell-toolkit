#ifndef DEBUG_HELPER
    #include <iostream>

    #define DEBUG_HELPER

    using namespace std;

    template<class T>
    void printArr(T *p, int n)
    {
        if(n==0) cout << "[]\n";
        cout << "[" << p[0];
        for(int i=1; i<n; ++i) cout << ", " << p[i];
        cout << "]" << endl;
    }

    template<class T>
    void printArr(T *p, int n, int m)
    {
        if(n==0) cout << "[]\n";
        cout << "[\n  " << p[0];
        for(int i=1; i<n*m; ++i) 
        {
            if( i%n == 0 ) cout << "\n";
            cout << ", " << p[i];
        }
        cout << "\n]" << endl;
    }
#endif