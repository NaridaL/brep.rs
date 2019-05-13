use std::fmt;
use std::marker;
use std::ops;
use num_traits::real::Real;
fn main() {
    println!("Hello, world!");

    let mut foos = V4::new(1.0, 2.0, 3.0, 4.0);
    foos += V4::new(0.0, 0.2, 0.4, 0.6);
    // let bar = 3.0 * V4::new(0.0, 0.2, 0.4, 0.6);
    // let car = bar * 3.0;
    // let dar = schur(bar, car);

    let a: f32 = 2.0;
    let b: f64 = 3.0;
    // let x = a * b;

    let x: V4<Foo> = V4::new(Foo {}, Foo {}, Foo {}, Foo {});
    let y = x;
    let z = x;

    let j = NVar("a".to_string());
    let l = Box::new(j) + Box::new(2.0);
}

trait Derivable: fmt::Display {
    fn derive(self, v: NVar) -> Box<Derivable>;
}
impl ops::Add<Box<Derivable>> for Box<Derivable> {
    type Output = Box<Addition>;
    fn add(self, rhs: Box<Derivable>) -> Box<Addition> {
        Box::new(Addition(self, rhs))
    }
}
impl ops::Add<Box<f64>> for Box<Derivable> {
    type Output = Box<Addition>;
    fn add(self, rhs: Box<f64>) -> Box<Addition> {
        Box::new(Addition(self, rhs))
    }
}

impl Derivable for f64 {
    fn derive(self, v: NVar) -> Box<Derivable> {
        Box::new(0.0)
    }
}

struct Addition(Box<Derivable>, Box<Derivable>);
impl fmt::Display for Addition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({} + {})", self.0, self.1)
    }
}
struct Multiplication(Box<Derivable>, Box<Derivable>);
impl fmt::Display for Multiplication {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({} * {})", self.0, self.1)
    }
}
#[derive(Debug, PartialEq)]
struct NVar(String);
impl fmt::Display for NVar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Derivable for NVar {
    fn derive(self, v: NVar) -> Box<Derivable> {
        if v == self {
            Box::new(1.0)
        } else {
            Box::new(0.0)
        }
    }
}

impl Derivable for Addition {
    fn derive(self, v: NVar) -> Box<Derivable> {
        Box::new(Addition(self.0.derive(v), self.1.derive(v)))
    }
}

#[derive(Copy, Clone)]
pub struct Foo {}

// struct V3<T>(T,  T, T);

// impl<T: ops::Div<Output = T> + std::marker::Copy> std::convert::From<V4<T>> for V3<T> {
//     fn from(error: V4<T>) -> Self {
//         V3(error.x / error.w, error.y / error.w, error.z / error.w)
//     }
// }
// pub fn schur<T: ops::Mul<Output = T>>(a: V4<T>, b: V4<T>) -> V4<T> {
//     V4::new(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w)
// }
// pub fn dot<T: ops::Mul<Output = S>, S: ops::Add<Output = S>>(a: V4<T>, b: V4<T>) -> S {
//     a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w
// }
// macro_rules! mul_primitive_V4_impl {
//     ($($t:ty)+) => ($(
//         impl ops::Mul<V4<$t>> for $t {
//             type Output = V4<$t>;
//                 fn mul(self, _rhs: V4<$t>) -> V4<$t> {
//                 V4 {
//                     x: self * _rhs.x,
//                     y: self * _rhs.y,
//                     z: self * _rhs.z,
//                     w: self * _rhs.w,
//                 }
//             }
//         }
//     )+)
// }
// mul_primitive_V4_impl! {usize u8 u16 u32 u64 u128 isize i8 i16 i32 i64 i128 f32 f64}

struct NURBS<T> {
    points: Vec<V4<T>>,
    knots: Vec<T>,
    degree: u32,
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_v4_add_assign() {
        let mut x = V4(1, 2, 3, 4);
        x += V4(1, 2, 3, 4);

        assert_eq!(x, V4(2, 4, 6, 8));
    }

    // #[test]
    // fn test_bad_add() {
    //     // This assert would fire and test will fail.
    //     // Please note, that private functions can be tested too!
    //     assert_eq!(bad_add(1, 2), 3);
    // }
}
