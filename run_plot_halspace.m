

h = HalfspaceRepresentation.getSquare();

figure(1);
clf();
xlim(2.0*[-1, 1]);
ylim(2.0*[-1, 1]);
axis square;
hold on;
h.plot()