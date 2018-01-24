(ns vmaths.core-test
  (:require [clojure.test :refer :all]
            [vmaths.core :refer :all]
            [clojure.core.matrix :as ma]))

(ma/set-current-implementation :vectorz)

(deftest vectors-test
  (let [v1 [4 15 16 18 32 0 39 32 15 19]
        v2 [7 38 29 23 11 32 25 5 34 17]]
    (testing "Vector"
      (is (= (mean v1) 19.0))
      (is (= (variance v1) 151.77777777777777))
      (is (= (pop-variance v1) 136.6))
      (is (= (ssd v1) 12.31981240838422))
      (is (= (psd v1) 11.687600266949584))
      (is (= (mad v1) 7.5))
      (is (= (meanad v1) 9.0))
      (mapv #(is (= (quantile v1 0.75 %1) %2))
            (range 1 10)
            [32.0 32.0 32.0 25.5 32.0 32.0 28.75 32.0 32.0])
      (is (= (median v2) 24.0))
      (is (= (mode v1) '(15 32)))
      (is (= (vrange v1) 39))
      (is (= (midrange v1) 19.5))
      (is (= (euclidean v1 v2) 59.0508255657785))
      (is (= (euclidean-squared v1 v2) 3487))
      (is (= (manhattan v1 v2) 159.0))
      (is (= (mae v1 v2) 15.9))
      (is (= (canberra v1 v2) 3.9976921256286575))
      (is (= (pearson v1 v2) -0.31231702416919843))
      (is (= (cosine v1 v2) 0.6888336157182211)))))

(deftest core-matrix-test
  (let [v1  (ma/array [4 15 16 18 32 0 39 32 15 19])
        v2 (ma/array [7 38 29 23 11 32 25 5 34 17])]
    (testing "Core matrix"
      (is (= (mean v1) 19.0))
      (is (= (variance v1) 151.77777777777777))
      (is (= (pop-variance v1) 136.6))
      (is (= (ssd v1) 12.31981240838422))
      (is (= (psd v1) 11.687600266949584))
      (is (= (mad v1) 7.5))
      (is (= (meanad v1) 9.0))
      (mapv #(is (= (quantile v1 0.75 %1) %2))
            (range 1 10)
            [32.0 32.0 32.0 25.5 32.0 32.0 28.75 32.0 32.0])
      (is (= (median v2) 24.0))
      (is (= (mode v1) [15.0 32.0]))
      (is (= (vrange v1) 39.0))
      (is (= (midrange v1) 19.5))
      (is (= (euclidean v1 v2) 59.0508255657785))
      (is (= (euclidean-squared v1 v2) 3487.0))
      (is (= (manhattan v1 v2) 159.0))
      (is (= (mae v1 v2) 15.9))
      (is (= (canberra v1 v2) 3.9976921256286575))
      (is (= (pearson v1 v2) -0.3123170241691985))
      (is (= (cosine v1 v2) 0.6888336157182211)))))

