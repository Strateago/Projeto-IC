from MatrixSpectralAnalysis import MatrixSpectralAnalysis

sa = MatrixSpectralAnalysis('Projeto-IC/arquivos/Problema_homogeneo_buck_40x40.mat')
method = 'MultiscaleILUfac'
av_error, res = sa.Solve(method)
sa.PlotSpectre(method, av_error)
sa.PlotResidues(method, res)