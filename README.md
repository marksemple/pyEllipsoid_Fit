# pyEllipsoid_Fit

A Python adaptation of Yury Petrov's [ellipsoid_fit](https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit) function for MATLAB. I use this for calibrating triaxis Magnetometers.

> Fits an ellipsoid or other conic surface into a 3D set of points approximating such a surface, allows some constraints, like orientation constraint and equal radii constraint. E.g., you can use it to fit a rugby ball, or a sphere. Returns both the algebraic description of the ellipsoid (the nine coefficients of the quadratic form) and the geometric description (center, radii, principal axes).

#### Notes to self for future updates
- Add example use cases (incl. sample data)
- Add unit testing
- Add a higher-level wrapper for better integration with Magnetometers
