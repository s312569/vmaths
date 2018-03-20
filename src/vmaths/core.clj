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

(defprotocol INormalisation
  (min-max [v])
  (standard-score [v])
  (positional-standardization [v])
  (unitization [v])
  (positional-unitization [v])
  (unitization-zero-min [v])
  (norm-in-range [v])
  (positional-norm-in-range [v])
  (normalization [v])
  (positional-norm [v])
  (norm-central-zero [v])
  (quotient-trans [v])
  (positional-quotient-trans [v])
  (quotient-trans-range [v])
  (quotient-trans-max [v])
  (quotient-trans-mean [v])
  (positional-quotient-trans-median [v])
  (quotient-trans-sum [v])
  (unit [v]))

(defprotocol INormalisable
  (do-normalise [o v rc])
  (normalisation-values [o f rc]))

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

(defn- mad-denominator
  [v]
  (let [m (mad v)]
    (if-not (= 0.0 m)
      (* 1.4826 m)
      (* 1.253314 (meanad v)))))

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
          (if (> (count v) 1)
            (let [{s :s m :m k :k} (vector-var v)]
              (-> (/ s (- k 1)) ma/sqrt))
            (println (str "WARNING: No ssd for vector of size " (count v)))))

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; normalise
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def default-vector-norm
  {:min-max
   (fn [v] (let [mi (ma/emin v) ma (ma/emax v)]
            [mi (- ma mi)]))
   
   :standard-score
   (fn [v] [(mean v) (ssd v)])

   :positional-standardization
   (fn [v] [(median v) (mad-denominator v)])

   :unitization
   (fn [v] [(mean v) (vrange v)])

   :positional-unitization
   (fn [v] [(median v) (vrange v)])

   :unitization-zero-min
   (fn [v] [(ma/emin v) (vrange v)])

   :norm-in-range
   (fn [v] (let [m (mean v)] [m (-> (ma/sub v m) ma/abs ma/emax)]))

   :positional-norm-in-range
   (fn [v] (let [m (median v)] [m (-> (ma/sub v m) ma/abs ma/emax)]))

   :quotient-trans
   (fn [v] [(ssd v)])

   :positional-quotient-trans
   (fn [v] [(mad-denominator v)])

   :quotient-trans-range
   (fn [v] [(vrange v)])

   :quotient-trans-max
   (fn [v] [(ma/emax v)])

   :quotient-trans-mean
   (fn [v] [(mean v)])

   :positional-quotient-trans-median
   (fn [v] [(median v)])

   :quotient-trans-sum
   (fn [v] [(ma/esum v)])

   :normalization
   (fn [v] (let [m (mean v)] [m (-> (ma/sub v m) ma/square ma/esum Math/sqrt)]))

   :positional-norm
   (fn [v] (let [m (median v)] [m (-> (ma/sub v m) ma/square ma/esum Math/sqrt)]))

   :norm-central-zero
   (fn [v] [(midrange v) (/ (vrange v) 2)])

   :unit
   (fn [v] [(mal/norm v)])})

(extend clojure.lang.PersistentVector INormalisation default-vector-norm)
(extend mikera.vectorz.Vector INormalisation default-vector-norm)
(extend mikera.vectorz.impl.StridedVector INormalisation default-vector-norm)
(extend mikera.vectorz.impl.ArraySubVector INormalisation default-vector-norm)

(defn- vector-norm
  [v [n d]]
  (if d
    (ma/emap #(-> (/ (- % n) d) double) v)
    (ma/emap #(-> (/ % n) double) v)))

(defn- norm-vals
  [v f]
  (let [[n d :as r] (f v)]
    (if d
      (if-not (zero? d) r (throw (Exception. "Zero denominator in normalisation.")))
      (if-not (zero? n) r (throw (Exception. "Zero denominator in normalisation."))))))

(def core-inormalisable
  {:do-normalise (fn [v [n d] _] (vector-norm v [n d]))
   :normalisation-values (fn [v f _] (norm-vals v f))})

(extend mikera.vectorz.Vector INormalisable core-inormalisable)
(extend mikera.vectorz.impl.StridedVector INormalisable core-inormalisable)
(extend mikera.vectorz.impl.ArraySubVector INormalisable core-inormalisable)

(extend-protocol INormalisable

  clojure.lang.PersistentVector

  (do-normalise [v cvs rc]
    (if (= (class (first v)) clojure.lang.PersistentVector)
      (condp = rc
        :rows (mapv #(do-normalise %1 %2 rc) v cvs)
        :cols (apply mapv vector (mapv #(do-normalise %1 %2 rc)
                                       (apply mapv vector v)
                                       cvs)))
      (vector-norm v cvs)))

  (normalisation-values [v f rc]
    (if (= (class (first v)) clojure.lang.PersistentVector)
      (if (= (-> (mapv count v) set count) 1)
        (mapv #(normalisation-values % f rc)
              (condp = rc
                :rows v
                :cols (apply mapv vector v)))
        (throw (Exception. "Unequal vector lengths")))
      (norm-vals v f))))

(extend-protocol INormalisable

  mikera.matrixx.Matrix

  (do-normalise [m cvs rc]
    (condp = rc
      :rows (ma/array (map #(do-normalise %1 %2 rc) (ma/rows m) cvs))
      :cols (->> (ma/array (map #(do-normalise %1 %2 rc) (ma/columns m) cvs))
                 ma/transpose)))

  (normalisation-values [m f rc]
    (condp = rc
      :cols (mapv #(normalisation-values % f rc) (ma/columns m))
      :rows (mapv #(normalisation-values % f rc) (ma/rows m)))))

(defn normalise
  "Normalises an object using function (nf) which can be one of:
  
  positional-unitization - ((x-median)/range)
  unit - (x/norm)
  positional-norm-in-range - ((x-median)/max(abs(x-median)))
  quotient-trans-sum - (x/sum)
  quotient-trans - (x/sd)
  quotient-trans-max - (x/max)
  positional-norm - ((x-median)/sqrt(sum((x-median)^2)))
  unitization - ((x-mean)/range)
  positional-standardization - ((x-median)/mad)
  norm-central-zero - ((x-midrange)/(range/2))
  quotient-trans-mean - (x/mean)
  quotient-trans-range - (x/range)
  positional-quotient-trans-median - (x/median)
  standard-score - ((x-mean)/sd)
  normalization - ((x-midrange)/(range/2))
  norm-in-range - ((x-mean)/max(abs(x-mean)))
  min-max - (x-minimum)/(maximum-minimum)
  unitization-zero-min - ((x-min)/range)
  positional-quotient-trans - (x/mad)

  Returns a vector containing normalised object and a vector of values
  used to normalize the object. For matrices the keyword 'dir'
  argument can be one of :rows or :cols; for vectors this argument is
  ignored."
  [o normfn & {:keys [dir] :or {dir :cols}}]
  (let [cvs (normalisation-values o normfn dir)]
    [(do-normalise o cvs dir) cvs]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; lazy seqs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- stream [f coll]
  (->> (if (:remaining coll) coll {:remaining coll})
       (iterate f)
       (take-while #(not (contains? % :end)))
       (map :yield)
       (filter #(not (nil? %)))))

(defn- stream-var
  [{[n & remaining] :remaining
    :keys [m s k]
    :or {m n s 0 k 0}}]
  (if-not n
    {:end true}
    (let [nk (inc k)
          nm (+ m (/ (- n m) nk))
          ns (+ s (* (- n m) (- n nm)))]
      {:m nm :k nk :s ns :yield {:m nm :s ns :k nk} :remaining remaining})))

(extend-protocol IStats

  clojure.lang.LazySeq

  (mean [l]
    (let [meanfn (fn [{[n & remaining] :remaining
                      :keys [i yield]
                      :or {i 0 yield 0}}]
                   (if-not n
                     {:end true}
                     (let [nm (/ (+ (* yield i) n) (+ i 1))]
                       {:remaining remaining :i (+ i 1) :yield nm})))]
      (map double (stream meanfn l))))

  (variance [l]
    (map (fn [{s :s m :m k :k}] (double (/ s (- k 1.0)))) (stream stream-var l)))

  (pop-variance [l]
    (map (fn [{s :s m :m k :k}] (double (/ s k))) (stream stream-var l)))

  (ssd [l]
    (map (fn [{s :s m :m k :k}] (double (-> (/ s (- k 1.0)) Math/sqrt)))
         (rest (stream stream-var l))))

  (psd [l]
    (map (fn [{s :s m :m k :k}] (double (-> (/ s k) Math/sqrt)))
         (rest (stream stream-var l))))

  (quantile [l p]
    (let [i (zipmap (range 1 6)
                    (mapv vector
                          (->> (take 5 l) sort vec)
                          [1.0 (+ 1 (* 2 p)) (+ 1 (* 4 p)) (+ 3 (* 2 p)) 5.0]
                          [0.0 (/ p 2) p (/ (+ 1 p) 2) 1.0]))
          f (fn [{[n & remaining] :remaining
                 :keys [s]}]
              (if n
                (let [[kn ns] (let [[fx f2 f3] (->> (vals s) first)
                                    [lx l2 l3] (->> (vals s) last)]
                                (cond (< n fx)
                                      [1 (assoc s 1 [n f2 f3])]
                                      (> n lx)
                                      [4 (assoc s 5 [n l2 l3])]
                                      :else
                                      [(->> (filter (fn [[[k [i x y]]
                                                         [k2 [i2 x2 y2]]]]
                                                      (and (>= n i) (< n i2)))
                                                    (partition-all 2 1 s))
                                            first first first)
                                       s]))
                      ns (->> (map (fn [[k v]] (if (> k kn) [(+ 1 k) v] [k v])) ns)
                              (map (fn [[k [i x y]]] [k [i (+ x y) y]]))
                              
                              (into {}))]
                  {:remaining remaining :s ns :yield [kn ns]})
                {:end true}))]
      (stream f {:remaining (drop 5 l) :s i}))))

(def qtest (lazy-seq [0.02 0.5 0.74 3.39 0.83 22.37 10.15 15.43 38.62 15.92
                      34.60 10.28 1.47 0.40 0.05 11.39 0.27 0.42 0.09 11.37]))
