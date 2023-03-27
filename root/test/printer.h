struct Printer {
    TCanvas canvas;
    std::string fname;
    int count;
    Printer(std::string fn)
      : canvas("test_raytiling", "Ray Tiling", 500, 500)
      , fname(fn)
      , count(0)
    {
        canvas.Print((fname + ".pdf[").c_str(), "pdf");
    }
    ~Printer() { canvas.Print((fname + ".pdf]").c_str(), "pdf"); }
    void operator()()
    {
        canvas.Print((fname + ".pdf").c_str(), "pdf");
        canvas.Print(Form("%s-%02d.png", fname.c_str(), count), "png");
        canvas.Print(Form("%s-%02d.svg", fname.c_str(), count), "svg");
        ++count;
    }
};

