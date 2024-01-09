@default_files = ('thesis.tex');

$pdf_mode = 1;

ensure_path('TEXINPUTS', './utils//',
            './chapters/introduction/figures',
            './chapters/polarons/figures',
            );
