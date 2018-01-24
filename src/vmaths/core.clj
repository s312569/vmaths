(ns vmaths.core
  (:require [clojure.core.matrix :as ma]
            [clojure.core.matrix.linear :as mal]
            [clojure.core.matrix.operators :as mos])
  (:import [org.apache.commons.math3.util Precision]))

(ma/set-current-implementation :vectorz)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; protocols
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defprotocol IStats
  (mean [v] "Returns the mean.")
  (variance [v] "Returns the variance.")
  (pop-variance [v] "Returns the population variance.")
  (ssd [v] "Sample standard deviation.")
  (psd [v] "Population standard deviation.")
  (mad [v] [v c] "Calculates the median absolute deviation about the
  median. Does not apply a distribution dependent scaling factor.")
  (meanad [v] "Calculates the mean absolute deviation about the
  median. Does not apply a distribution dependent scaling factor.")
  (quantile [v p] [v p m] "Returns the pth quantile")
  (median [v] "Returns the median.")
  (mode [v] "Returns the mode.")
  (vrange [v] "Returns the range of a vector.")
  (midrange [v] "Returns the mid-range of a vector."))

(defprotocol IDistance
  (euclidean [a b] "Returns the euclidean distance.")
  (euclidean-squared [a b] "Returns the euclidean squared distance.")
  (manhattan [a b] "Returns the manhattan distance.")
  (mae [a b] "Returns the mae distance.")
  (canberra [a b] "Returns the canberra distance."))

(defprotocol ISimilarity
  (pearson [a b])
  (cosine [a b]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn float=
  ([^double x ^double y ^double nepsilons]
   (Precision/equals x y (* Precision/EPSILON nepsilons)))
  ([^double x ^double y]
   (Precision/equals x y Precision/EPSILON)))

(defn float<
  ([^double x ^double y ^double nepsilons]
   (< (Precision/compareTo x y (* Precision/EPSILON nepsilons)) 0))
  ([^double x ^double y]
   (< (Precision/compareTo x y Precision/EPSILON) 0)))

(defn float<=
  ([^double x ^double y ^double nepsilons]
   (<= (Precision/compareTo x y (* Precision/EPSILON nepsilons)) 0))
  ([^double x ^double y]
   (<= (Precision/compareTo x y Precision/EPSILON) 0)))

(defn float>
  ([^double x ^double y ^double nepsilons]
   (> (Precision/compareTo x y (* Precision/EPSILON nepsilons)) 0))
  ([^double x ^double y]
   (> (Precision/compareTo x y Precision/EPSILON) 0)))

(defn float>=
  ([^double x ^double y ^double nepsilons]
   (>= (Precision/compareTo x y (* Precision/EPSILON nepsilons)) 0))
  ([^double x ^double y]
   (>= (Precision/compareTo x y Precision/EPSILON) 0)))

(defn- ma-count
  [v]
  (-> (ma/shape v) first))

(defn- ma-last
  [v]
  (ma/mget v (- (ma-count v) 1)))

(defn- vector-var
  [v]
  (reduce (fn [{m :m s :s k :k} x]
            (let [nk (inc k)
                  nm (+ m (/ (- x m) nk))]
              {:m nm :k nk :s (+ s (* (- x m) (- x nm)))}))
          {:m (ma/mget v 0) :s 0 :k 0}
          v))

(defn- calc-quant-dis
  [p v m]
  (if (#{mikera.vectorz.Vector mikera.vectorz.impl.StridedVector
         mikera.vectorz.impl.ArraySubVector clojure.lang.PersistentVector}
       (class v))
    [::default m]
    [(class v) m]))

(defmulti calc-quantile #'calc-quant-dis)

(defmethod calc-quantile [::default 1]
  [p v _]
  (let [np (* p (ma-count v))
        j (int (ma/floor np))
        ij (- j 1)]
    (if (float= (- np j) 0.0)
      (ma/mget v ij)
      (ma/mget v (+ ij 1)))))

(defmethod calc-quantile [::default 2]
  [p v _]
  (cond (float= 0.0 p)
        (ma/mget v 0)
        (float= 1.0 p)
        (ma-last v)
        :else
        (let [np (* p (ma-count v))
              j (int (ma/floor np))
              ij (- j 1)]
          (if (float= (- np j) 0)
            (* (+ (ma/mget v ij) (ma/mget v (+ 1 ij))) 0.5)
            (ma/mget v (+ ij 1))))))

(defmethod calc-quantile [::default 3]
  [p v _]
  (let [c (ma-count v)]
    (if (<= p (/ 0.5 c))
      (ma/mget v 0)
      (let [np (* p c)
            j (int (Math/floor np))
            ij (if (float= (- np j) 0.50)
                 (if (= 0 (mod j 2)) (- j 1) j)
                 (- (Math/round np) 1))]
        (ma/mget v ij)))))

(defn- quant-helper
  [np v]
  (let [j (Math/floor np)
        ij (- j 1)]
    (+ (ma/mget v ij) (* (- np j) (- (ma/mget v j) (ma/mget v ij))))))

(defmethod calc-quantile [::default 4]
  [p v _]
  (let [c (ma-count v)]
    (cond (float< p (/ 1 c)) (ma/mget v 0)
          (float= 1.0 p)(ma-last v)
          :else (quant-helper (* p c) v))))

(defmethod calc-quantile [::default 5]
  [p v _]
  (let [c (ma-count v)]
    (cond (float< p (/ 0.5 c)) (ma/mget v 0)
          (float>= p (/ (- c 0.5) c)) (ma-last v)
          :else (quant-helper (+ (* p c) 0.5) v))))

(defmethod calc-quantile [::default 6]
  [p v _]
  (let [c (ma-count v)]
    (cond (float< p (/ 1 (+ c 1))) (ma/mget v 0)
          (float>= p (/ c (+ 1 c))) (ma-last v)
          :else (quant-helper (* p (+ 1 c)) v))))

(defmethod calc-quantile [::default 7]
  [p v _]
  (let [c (ma-count v)]
    (cond (float= p 1.0) (ma-last v)
          :else (quant-helper (+ (* (- c 1) p) 1) v))))

(defmethod calc-quantile [::default 8]
  [p v _]
  (let [c (ma-count v)
        t1 (/ 1 3)
        t2 (/ 2 3)]
    (cond (float< p (/ t2 (+ c t1))) (ma/mget v 0)
          (float>= p (/ (- c t1) (+ c t2))) (ma-last v)
          :else (quant-helper (+ (* (+ c t1) p) t1) v))))

(defmethod calc-quantile [::default 9]
  [p v _]
  (let [c (ma-count v)
        t5 (/ 5 8)
        t4 (/ 1 4)
        t3 (/ 3 8)]
    (cond (float< p (/ t5 (+ c t4))) (ma/mget v 0)
          (float>= p (/ (- c t3) (+ c t4))) (ma-last v)
          :else (quant-helper (+ (* (+ c t4) p) t3) v))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; vectors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def default-stats
  {:mean (fn [v]
           (double (/ (ma/ereduce + v) (ma-count v))))

   :variance (fn [v]
               (let [{s :s m :m k :k} (vector-var v)]
                 (double (/ s (- k 1.0)))))

   :pop-variance (fn [v]
                   (let [{s :s m :m k :k} (vector-var v)]
                     (double (/ s k))))

   :ssd (fn [v]
          (let [{s :s m :m k :k} (vector-var v)]
            (-> (/ s (- k 1)) ma/sqrt)))

   :psd (fn [v]
          (let [{s :s m :m k :k} (vector-var v)]
            (-> (/ s k) ma/sqrt)))

   :mad (fn [v]
          (let [m (median v)] (-> (ma/emap #(ma/abs (- % m)) v) median)))

   :meanad (fn [v]
             (let [m (median v)] (/ (->> (ma/emap #(ma/abs (- % m)) v)
                                         (ma/ereduce +))
                                    (ma-count v))))

   :quantile (fn ([v p] (double (calc-quantile p (-> (sort < v) vec) 7)))
               ([v p m] (double (calc-quantile p (-> (sort < v) vec) m))))
   
   :median (fn [v]
             (double (calc-quantile 0.5 (-> (sort < v) vec) 2)))

   :mode (fn [v]
           (->> (frequencies v) (group-by val) (sort-by key) last val (mapv key)))

   :vrange (fn [v]
             (- (ma/emax v) (ma/emin v)))

   :midrange (fn [v]
               (/ (+ (ma/emax v) (ma/emin v)) 2.0))})

(extend clojure.lang.PersistentVector IStats default-stats)
(extend mikera.vectorz.Vector IStats default-stats)
(extend mikera.vectorz.impl.StridedVector IStats default-stats)
(extend mikera.vectorz.impl.ArraySubVector IStats default-stats)

(def default-distance
  {:euclidean (fn [a b]
                (ma/distance a b))

   :euclidean-squared (fn [a b]
                        (->> (ma/sub a b) ma/square ma/esum))

   :manhattan (fn [a b]
                (->> (ma/sub a b) ma/abs ma/esum))

   :mae (fn [a b]
          (double (/ (->> (ma/sub a b) ma/abs ma/esum) (ma-count a))))

   :canberra (fn [a b]
               (->> (ma/emap (fn [x y]
                               (if (and (float= x 0.0) (float= y 0.0))
                                 0.0
                                 (double
                                  (/ (-> (- x y) ma/abs) (+ (ma/abs x) (ma/abs y))))))
                             a b)
                    ma/esum))})

(extend clojure.lang.PersistentVector IDistance default-distance)
(extend mikera.vectorz.Vector IDistance default-distance)
(extend mikera.vectorz.impl.StridedVector IDistance default-distance)
(extend mikera.vectorz.impl.ArraySubVector IDistance default-distance)

(def default-sim
  {:pearson (fn [a b]
              (let [c (ma-count a)
                    sumx (ma/esum a)
                    sumy (ma/esum b)
                    n (- (-> (ma/mmul a b) ma/mget) (/ (* sumx sumy) c))
                    d (* (ma/sqrt
                          (- (-> (ma/square a) ma/esum)
                             (/ (* sumx sumx) c)))
                         (ma/sqrt
                          (- (-> (ma/square b) ma/esum)
                             (/ (* sumy sumy) c))))]
                (if (> d 0.0) (double (/ n d)) 0.0)))

   :cosine (fn [a b]
             (let [c (ma-count a)
                   n (-> (ma/mmul a b) ma/mget)
                   d (* (ma/sqrt (-> (ma/square a) ma/esum))
                        (ma/sqrt (-> (ma/square b) ma/esum)))]
               (if (> d 0.0) (double (/ n d)) 0.0)))})

(extend clojure.lang.PersistentVector ISimilarity default-sim)
(extend mikera.vectorz.Vector ISimilarity default-sim)
(extend mikera.vectorz.impl.StridedVector ISimilarity default-sim)
(extend mikera.vectorz.impl.ArraySubVector ISimilarity default-sim)


