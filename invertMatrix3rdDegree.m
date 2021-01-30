function Aodwr = invertMatrix3rdDegree(A)

detA = det(A);

Aodwr = [
    (A(2,2)*A(3,3)-A(2,3)*A(3,2))/detA (A(1,3)*A(3,2)-A(1,2)*A(3,3))/detA (A(1,2)*A(2,3)-A(1,3)*A(2,2))/detA
    (A(2,3)*A(3,1)-A(2,1)*A(3,3))/detA (A(1,1)*A(3,3)-A(1,3)*A(3,1))/detA (A(1,3)*A(2,1)-A(1,1)*A(2,3))/detA
    (A(2,1)*A(3,2)-A(2,2)*A(3,1))/detA (A(1,2)*A(3,1)-A(1,1)*A(3,2))/detA (A(1,1)*A(2,2)-A(1,2)*A(2,1))/detA
];

%Aodwr = inv(A); 