
;;; trie.el --- ternary search tree package


;; Copyright (C) 2004-2007 Toby Cubitt

;; Author: Toby Cubitt <toby-predictive@dr-qubit.org>
;; Version: 0.1
;; Keywords: trie, ternary search tree, completion
;; URL: http://www.dr-qubit.org/emacs.php


;; This file is NOT part of Emacs.
;;
;; This program is free software; you can redistribute it and/or
;; modify it under the terms of the GNU General Public License
;; as published by the Free Software Foundation; either version 2
;; of the License, or (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program; if not, write to the Free Software
;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
;; MA 02110-1301, USA.


;;; Commentary:
;;
;; A ternary search tree associates data with keys. The keys can be any
;; ordered sequence of elements: vector, list or string. It stores them
;; in such a way that both storage size and data lookup are reasonably
;; space- and time-efficient, respectively. But more importantly,
;; returning all strings with a given prefix in alphabetical or any other
;; sort-order is also time-efficient.
;;
;; Internally, a ternary search tree consists of two elements, the root
;; node and the comparison function. The latter must take two arguments
;; of the same type as individual elements of a key, and must return nil if
;; the first is smaller than the second.
;;
;; Note that there are two uses for a ternary search tree: as a
;; "dictionary", in which case only presence of absence of a key is
;; significant, or as an associative array, in which case keys are
;; associated with data. Other similar data types often only implement
;; the first of these, leaving it up to you to implement an associative
;; on top of this. However, this would be somewhat space-inefficient in a
;; ternary search tree, so this package does the opposite: it implements
;; associative arrays, and leaves it up to you to use them as
;; dictionaries if you so desire.
;;
;; You create a ternary search tree using `trie-create', create an
;; association using `trie-insert', retrieve an association using
;; `trie-member', find completions of a sequence using
;; `trie-complete', find completions and sort them in a specified order
;; using `trie-complete-ordered', and map over a tree using
;; `trie-map' or `trie-mapcar'.
;;
;;
;; This package uses the AVL tree package avl-tree.el and the ternary
;; heap package heap.el.


;;; Change Log:
;;
;; Version 0.1
;; * Initial release (complete rewrite from scratch of tstree.el!)
;; * Ternary search trees are now implemented as a tree of avl trees, which
;;   has numerous advantages: self-balancing trees guarantee O(log n)
;;   complexity regardless of how the tree is built; deletion is now done
;;   correctly.
;; * Up to "tstree"->"trie" renaming, functions are almost drop-in
;;   replacements for tstree.el functions. However, insertion and rank
;;   functions are no longer stored in the data structure, so corresponidng
;;   arguments are no longer optional. And functions can no longer operate
;;   over multiple data structures at once; i.e. they no longer accept lists
;;   of trees or prefixes as arguments. (These features belong in higher level
;;   libraries, and the efficiency loss is negligible.)
;; * `trie.el' is now general enough to implement all sorts of tries, not just
;;   ternary search trees (though these remain the default).



;;; Code:

(provide 'trie)
(require 'avl-tree)
(require 'heap)




;;; ================================================================
;;;                Replacements for CL functions

;; copied from cl-extra.el
(defun trie--subseq (seq start &optional end)
  "Return the subsequence of SEQ from START to END.
If END is omitted, it defaults to the length of the sequence.
If START or END is negative, it counts from the end."
  (if (stringp seq) (substring seq start end)
    (let (len)
      (and end (< end 0) (setq end (+ end (setq len (length seq)))))
      (when (< start 0)
	(setq start (+ start (or len (setq len (length seq))))))
      (cond ((listp seq)
	     (if (> start 0) (setq seq (nthcdr start seq)))
	     (if end
		 (let ((res nil))
		   (while (>= (setq end (1- end)) start)
		     (push (pop seq) res))
		   (nreverse res))
	       (copy-sequence seq)))
	    (t
	     (or end (setq end (or len (length seq))))
	     (let ((res (make-vector (max (- end start) 0) nil))
		   (i 0))
	       (while (< start end)
		 (aset res i (aref seq start))
		 (setq i (1+ i) start (1+ start)))
	       res))))))


(defun trie--seq-append (seq el)
  "Append EL to the end of sequence SEQ."
  (cond
   ((stringp seq) (concat seq (string el)))
   ((vectorp seq) (vconcat seq (vector el)))
   ((listp seq)	  (append seq (list el)))))



;;; ================================================================
;;;  Internal functions only for use within the trie package

;;; ----------------------------------------------------------------
;;; Functions and macros for handling a trie.

(defstruct
  (trie-
   :named
   (:constructor nil)
   (:constructor trie--create
		 (comparison-function &optional type
		  &aux
		  (cmpfun `(lambda (a b)
			     (setq a (trie--node-split a)
				   b (trie--node-split b))
			     (cond ((eq a 'trie--terminator)
				    (if (eq b 'trie--terminator) nil t))
				   ((eq b 'trie--terminator) nil)
				   (t (,comparison-function a b)))))
		  (createfun
		   (cond
		    ((or (eq type 'avl) (null type)) 'avl-tree-create)
		    (t (error "trie--create: unrecognized TYPE"))))
		  (insertfun
		   (cond
		    ((or (eq type 'avl) (null type)) 'avl-tree-enter)
		    (t (error "trie--create: unrecognized TYPE"))))
		  (deletefun
		   (cond
		    ((or (eq type 'avl) (null type)) 'avl-tree-delete)
		    (t (error "trie--create: unrecognized TYPE"))))
		  (retrievefun
		   (cond
		    ((or (eq type 'avl) (null type)) 'avl-tree-member)
		    (t (error "trie--create: unrecognized TYPE"))))
		  (mapfun
		   (cond
		    ((or (eq type 'avl) (null type)) 'avl-tree-mapc)
		    (t (error "trie--create: unrecognized TYPE"))))
		  (emptyfun
		   (cond
		    ((or (eq type 'avl) (null type)) 'avl-tree-empty)
		    (t (error "trie--create: unrecognized TYPE"))))
		  (root (trie--node-create-root createfun cmpfun))
		  ))
   (:copier nil))
  root comparison-function cmpfun
  createfun insertfun deletefun retrievefun mapfun emptyfun)


(defmacro trie--wrap-cmpfun (cmpfun)
  ;; wrap CMPFUN for use in a subtree
  `(lambda (a b)
     (setq a (trie--node-split a)
	   b (trie--node-split b))
     (cond ((eq a 'trie--terminator)
	    (if (eq b 'trie--terminator) nil t))
	   ((eq b 'trie--terminator) nil)
	   (t (,cmpfun a b)))))


;; (defmacro trie--wrap-cmpfun (cmpfun)
;;   ;; Wrap comparison function so that 'trie--terminator is always less than
;;   ;; anything other than another 'trie--terminator, and it unpacks the split
;;   ;; cell from a trie node.
;;   `(lambda (a b)
;;      (setq a (trie--node-split a)
;; 	   b (trie--node-split b))
;;      (cond ((eq a 'trie--terminator)
;; 	    (if (eq b 'trie--terminator) nil t))
;; 	   ((eq b 'trie--terminator) nil)
;; 	   (t (,cmpfun a b)))))


;;; ----------------------------------------------------------------
;;; Functions and macros for handling a trie node.

(defstruct
  (trie--node
   (:type vector)
   (:constructor nil)
   (:constructor trie--node-create
		 (split tree
		  &aux (subtree (funcall (trie--createfun tree)
					 (trie--cmpfun tree)))))
   (:constructor trie--node-create-data
		 (data &aux (split 'trie--terminator) (subtree data)))
   (:constructor trie--node-create-dummy
		 (split &aux (subtree nil)))
   (:constructor trie--node-create-root
		 (createfun cmpfun
		  &aux (split nil) (subtree (funcall createfun cmpfun))))
   (:copier nil))
   split subtree)

;; data is stored in the subtree cell of a terminal node
(defalias 'trie--node-data 'trie--node-subtree)

(defsetf trie--node-data (node)
  `(setf (trie--node-subtree ,node)))


(defun trie--node-find (tree sequence)
  ;; Returns the node corresponding to SEQUENCE, or nil if none found.
  (let ((node (trie--root tree))
	(len (length sequence))
	(i -1))
    ;; descend tree until we find SEQUENCE or run out of tree
    (while (and node (< (incf i) len))
      (setq node
	    (funcall (trie--retrievefun tree)
		     (trie--node-subtree node)
		     (trie--node-create-dummy (elt sequence i)))))
    node))


(defmacro trie--node-find-data (node)
  ;; Return data node from NODE's subtree, or nil if NODE has no data node in
  ;; its subtree.
  `(funcall (trie--retrievefun tree)
	    (trie--node-subtree ,node)
	    (trie--node-create-dummy 'trie--terminator)))


(defmacro trie--find-data (node)
  ;; Return data associated with sequence corresponding to NODE, or nil if
  ;; sequence has no associated data.
  `(let ((node (trie--node-find-data ,node)))
     (trie--node-data node)))


(defmacro trie--node-data-p (node)
  ;; Return t if NODE is a data node, nil otherwise.
  `(eq (trie--node-split ,node) 'trie--terminator))


(defun trie--mapc (trie--mapc--function trie--mapc--mapfun
					    trie--root seq
					    &optional reverse)
  ;; Apply FUNCTION to all elements in a trie beneath ROOT, which should
  ;; correspond to the sequence SEQ. TRIE-FUNCTION is passed two arguments:
  ;; the trie node itself and the sequence it corresponds to. It is applied
  ;; in ascending order, or descending order if REVERSE is non-nil. ARGS are
  ;; passed as additional arguments to TRIE-FUNCTION
  (funcall
   trie--mapc--mapfun
   (lambda (node)
     ;; data node: apply function
     (if (trie--node-data-p node)
	 (funcall trie--mapc--function node seq)
       ;; internal node: append split value to seq and keep descending
       (trie--mapc trie--mapc--function trie--mapc--mapfun node
		     (trie--seq-append (copy-sequence seq)
					 (trie--node-split node))
		     reverse)))
   ;; TRIE--TREE--MAPFUN target
   (trie--node-subtree trie--root)
   reverse))



;;; ================================================================
;;;    The public functions which operate on ternary search trees.

(defalias 'trie-create 'trie--create
  "Return a new ternary search tree that uses comparison function CMPFUN.

A ternary search tree stores sequences (string, vector or list)
along with associated data. CMPFUN should accept two arguements,
each an individual element of such a sequence, and return t if
the first is smaller than the second.")


(defalias 'trie-comparison-function 'trie--comparison-function
  "Return the comparison function for the ternary search tree TREE.")


(defalias 'trie-p 'trie--p
  "Return t if TREE is a ternary search tree, nil otherwise.")


(defun trie-empty (tree)
  "Return t if the ternary search tree TREE is empty, nil otherwise."
  (funcall (trie--emptyfun tree)
	   (trie--node-subtree (trie--root tree))))


(defun trie-map (function tree &optional type reverse)
  "Modify all elements in ternary search tree TREE
by applying FUNCTION to them.

FUNCTION should take two arguments: a sequence stored in the tree
and its associated data. Its return value replaces the existing
data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

FUNCTION is applied in ascending order, or descending order if
REVERSE is non-nil."
  (let ((trie-mapc--function function)) ; try to avoid dynamic binding bugs
    (trie--mapc
     (lambda (node seq)
       (setf (trie--node-data node)
	     (funcall trie-mapc--function seq (trie--node-data node))))
     (trie--mapfun tree)
     (trie--root tree)
     (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
     reverse)))


(defun trie-mapc (function tree &optional type reverse)
  "Apply FUNCTION to all elements in ternary search tree TREE for
side effect only.

FUNCTION should take two arguments: a sequence stored in the tree
and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

FUNCTION is applied in ascending order, or descending order if
REVERSE is non-nil."
  (let ((trie-mapc--function function)) ; try to avoid dynamic binding bugs
    (trie--mapc
     (lambda (node seq)
       (funcall trie-mapc--function seq
		(trie--node-data node)))
     (trie--mapfun tree)
     (trie--root tree)
     (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
     reverse)))


(defun trie-mapf (function combinator tree &optional type reverse)
  "Apply FUNCTION to all elements in ternary search tree TREE,
and combine the results using COMBINATOR.

FUNCTION should take two arguments: a sequence stored in the tree
and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

FUNCTION is applied in ascending order, or descending order if
REVERSE is non-nil."
  (let ((trie-mapc--function function) ; try to avoid dynamic binding bugs
	trie-mapcar--accumulate)
    (trie--mapc
     (lambda (node seq)
       (funcall combinator
		(funcall trie-mapc--function seq (trie--node-data node))
		trie-mapcar--accumulate))
     (trie--mapfun tree)
     (trie--root tree)
     (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
     reverse)))


(defun trie-mapcar (function tree &optional type reverse)
  "Apply FUNCTION to all elements in ternary search tree TREE,
and make a list of the results.

FUNCTION should take two arguments: a sequence stored in the tree
and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

FUNCTION is applied in ascending order, or descending order if
REVERSE is non-nil."
  (trie-mapf function 'cons tree type reverse))


;; ----------------------------------------------------------------
;;                        Inserting data

(defun trie-insert (tree key &optional data updatefun)
  "Associate DATA with KEY in ternary search tree TREE.

If KEY already exists in TREE, then DATA replaces the existing
association, unless UPDATEFUN is supplied. Note that if DATA is
*not* supplied, this means that the existing association of KEY
will be replaced by nil.

If UPDATEFUN is supplied and KEY already exists in TREE,
UPDATE-FUNCTION is called with two arguments: DATA and the
existing association of KEY. Its return value becomes the new
association for KEY.

Returns the new association of KEY."
  (let ((node (trie--root tree))
	(len (length key))
	(keep-old (lambda (a b) b))
	update-old
	(i -1))
    ;; Descend tree, adding nodes for non-existent elements of KEY. The update
    ;; function passed to `trie--enterfun' ensures that existing nodes are
    ;; left intact.
    (while (< (incf i) len)
      (setq node (funcall (trie--insertfun tree)
			  (trie--node-subtree node)
			  (trie--node-create (elt key i) tree)
			  keep-old)))
    ;; If UPDATEFUN was supplied, wrap it for passing to `trie--enterfun'.
    (if updatefun
	(setq update-old
	      (lambda (a b)
		(setf (trie--node-data b)
		      (funcall updatefun
			       (trie--node-data a)
			       (trie--node-data b))))))
    ;; Create or update data node.
    (setq node (funcall (trie--insertfun tree)
			(trie--node-subtree node)
			(trie--node-create-data data)
			update-old))
    (trie--node-data node)))  ; return new data



;; ----------------------------------------------------------------
;;                        Deleting data

(defun trie-delete (tree key)
  "Delete KEY and its associated data from TREE.

If KEY was deleted, a cons cell containing KEY and its
association is returned. Returns nil if KEY does not exist in
TREE."
  (let (trie--deleted-node)
    (declare (special trie--deleted-node))
    (trie--do-delete (trie--root tree) key
		       (trie--deletefun tree) (trie--emptyfun tree))
    (cons key (trie--node-data trie--deleted-node))))


(defun trie--do-delete (node seq deletefun emptyfun)
  ;; Delete SEQ starting from trie node NODE, and return non-nil if we
  ;; deleted a node.
  (declare (special trie--deleted-node))
  ;; if SEQ is empty, try to delete data node and return non-nil if we did
  ;; (return value of DELETEFUN is the deleted data, which is always non-nil
  ;; for a trie)
  (if (= (length seq) 0)
      (setq trie--deleted-node
	    (funcall deletefun
		     (trie--node-subtree node)
		     (trie--node-create-dummy 'trie--terminator)))
    ;; otherwise, delete on down (return value of DELETEFUN is the deleted
    ;; data, which is always non-nil for a trie)
    (funcall deletefun
	     (trie--node-subtree node)
	     (trie--node-create-dummy (elt seq 0))
	     (lambda (n)
	       (and (trie--do-delete n (trie--subseq seq 1)
				       deletefun emptyfun)
		    (funcall emptyfun (trie--node-subtree n)))))))


;; ----------------------------------------------------------------
;;                       Retrieving data

(defun trie-member (tree key &optional nilflag)
  "Return the data associated with KEY in the tree TREE,
or nil if KEY does not exist in TREE.

Optional argument NILFLAG specifies a value to return instead of
nil if KEY does not exist in TREE. This allows a non-existent KEY
to be distinguished from an element with a null association. (See
also `trie-member-p', which does this for you.)"
  ;; find node corresponding to key, then find data node, then return data
  (let (node)
    (or (and (setq node (trie--node-find tree key))
	     (trie--find-data node))
	nilflag)))

(defun trie-member-p (tree key)
  "Return t if KEY exists in TREE, nil otherwise."
  (let ((flag '(nil))) (not (eq flag (trie-member tree key flag)))))




;; ----------------------------------------------------------------
;;                         Completing

;; (defun trie--do-complete (node seq accumulator store filter)
;;   ;; Return all elements in a trie beneath NODE, which must correspond to
;;   ;; the sequence SEQ. More specifically, if an element passes the FILTER
;;   ;; function, ACCUMULATOR is called with two arguments: a cons of the
;;   ;; sequence corresponding to that element and its associated data, and
;;   ;; STORE. The elements are traversed in order.
;;   (avl-tree-mapc
;;    (lambda (node)
;;      ;; data node: found potential completion
;;      (if (trie--node-data-p node)
;; 	 (let ((data (trie--node-data node)))
;; 	   ;; if it passes filter, add it to completions list
;; 	   (if (or (null filter) (funcall filter seq data))
;; 	       (funcall accumulator (cons seq data) store)))
;;        ;; internal node: append split value to seq and keep descending
;;        (trie--do-complete
;; 	node
;; 	(cond
;; 	 ((stringp seq)
;; 	  (concat (copy-sequence seq)
;; 		  (string (trie--node-split node))))
;; 	 ((vectorp seq)
;; 	  (vconcat (copy-sequence seq)
;; 		   (vector (trie--node-split node))))
;; 	 ((listp seq)
;; 	  (append (copy-sequence seq)
;; 		  (list (trie--node-split node))))
;; 	 (t (error "trie-complete: invalid KEY type, sequencep")))
;; 	 accumulator store filter)))
;;    ;; avl-tree-mapc target
;;    (trie--node-subtree node)))


;; (defun trie-complete (tree key &optional maxnum filter)
;;   "Return an alist containing all completions of KEY
;; in ternary searh tree TREE along with their associated data, in
;; \"lexical\" order (i.e. the order defined by the tree's
;; comparison function). If no completions are found, return nil.

;; KEY must be a sequence (vector, list or string) containing
;; elements of the type used to reference data in the tree. (If KEY
;; is a string, it must be possible to apply `string' to individual
;; elements of the sequences stored in the tree.) The completions
;; returned in the alist will be sequences of the same type as KEY.

;; The optional integer argument MAXNUM limits the results to the
;; first MAXNUM completions. Otherwise, all completions are
;; returned.

;; The FILTER argument sets a filter function for the
;; completions. If supplied, it is called for each possible
;; completion with two arguments: the completion, and its associated
;; data. If the filter function returns nil, the completion is not
;; included in the results, and does not count towards MAXNUM."

;;   (let* ((node (trie--node-find tree key))
;; 	 (trie--stack (make-vector 1 nil))
;; 	 (accumulator
;; 	  (lambda (a stk)
;; 	    (message "%s" a)
;; 	    (message "%s" stk)
;; 	    (aset stk 0 (cons a (aref stk 0)))
;; 	    (and maxnum
;; 		 (>= (length (aref stk 0)) maxnum)
;; 		 (throw 'trie-complete--done nil)))))
;;     (when node
;;       (catch 'trie-complete--done
;; 	(trie--do-complete node key accumulator trie--stack filter))
;;       (nreverse (aref trie--stack 0)))))


;; (defun trie-complete-ordered (tree key rankfun &optional maxnum filter)
;;   "Return an alist containing all completions of SEQUENCE found
;; in ternary searh tree TREE along with their associated data, in
;; the order defined by RANKFUN. If no completions are found, return
;; nil.

;; KEY must be a sequence (vector, list or string) containing
;; elements of the type used to reference data in the tree. (If KEY
;; is a string, it must be possible to apply `string' to individual
;; elements of the sequences stored in the tree.) The completions
;; returned in the alist will be sequences of the same type as KEY.

;; RANKFUN must accept two arguments, both cons cells. The car
;; contains a sequence from the tree (of the same type as KEY), the
;; cdr contains its associated data.

;; The optional integer argument MAXNUM limits the results to the
;; first MAXNUM completions. Otherwise, all completions are
;; returned.

;; The FILTER argument sets a filter function for the
;; completions. If supplied, it is called for each possible
;; completion with two arguments: the completion, and its associated
;; data. If the filter function returns nil, the completion is not
;; included in the results, and does not count towards MAXNUM."

;;   (let ((node (trie--node-find tree key))
;; 	(heap (heap-create
;; 	       `(lambda (a b) (not (,rankfun a b)))
;; 	       (when maxnum (1+ maxnum))))
;; 	(accumulator (lambda (a hp)
;; 		       (heap-add hp a)
;; 		       (when (and maxnum (> (heap-size hp) maxnum))
;; 			 (heap-delete-root hp)))))
;;     (when node (trie--do-complete node key accumulator heap filter))

;;     (let (completions)
;;       (while (not (heap-empty heap))
;; 	(setq completions (cons (heap-delete-root heap) completions)))
;;       completions)))


(defun trie-complete (tree key &optional maxnum filter)
  "Return an alist containing all completions of KEY
in ternary searh tree TREE along with their associated data, in
\"lexical\" order (i.e. the order defined by the tree's
comparison function). If no completions are found, return nil.

KEY must be a sequence (vector, list or string) containing
elements of the type used to reference data in the tree. (If KEY
is a string, it must be possible to apply `string' to individual
elements of the sequences stored in the tree.) The completions
returned in the alist will be sequences of the same type as KEY.

The optional integer argument MAXNUM limits the results to the
first MAXNUM completions. Otherwise, all completions are
returned.

The FILTER argument sets a filter function for the
completions. If supplied, it is called for each possible
completion with two arguments: the completion, and its associated
data. If the filter function returns nil, the completion is not
included in the results, and does not count towards MAXNUM."

  (let* ((node (trie--node-find tree key))
	 (trie--stack (make-vector 1 nil))
	 (accumulator
	  (lambda (node seq)
	    (let ((data (trie--node-data node)))
	      (when (or (null filter) (funcall filter seq data))
		(aset trie--stack 0
		      (cons (cons seq data) (aref trie--stack 0)))
		(and maxnum
		     (>= (length (aref trie--stack 0)) maxnum)
		     (throw 'trie-complete--done nil)))))))
    (when node
      (catch 'trie-complete--done
	(trie--mapc accumulator (trie--mapfun tree) node key t))
      (aref trie--stack 0))))


;; Note: We use a partial heap-sort to find the k=MAXNUM highest ranked
;; completions among n possibilities. This has worst-case time complexity
;; O(n log k), and is both simple and elegant. An optimal algorithm
;; (e.g. partial quick-sort where the irrelevant partition is discarded
;; at each step) would have complexity O(n + k log k), but is probably
;; not worth the extra coding effort, and would have worse space
;; complexity. (I haven't done any benchmarking, though, so feel free to
;; do so and let me know the results!)
(defun trie-complete-ordered (tree key rankfun &optional maxnum filter)
  "Return an alist containing all completions of SEQUENCE found
in ternary searh tree TREE along with their associated data, in
the order defined by RANKFUN. If no completions are found, return
nil.

KEY must be a sequence (vector, list or string) containing
elements of the type used to reference data in the tree. (If KEY
is a string, it must be possible to apply `string' to individual
elements of the sequences stored in the tree.) The completions
returned in the alist will be sequences of the same type as KEY.

RANKFUN must accept two arguments, both cons cells. The car
contains a sequence from the tree (of the same type as KEY), the
cdr contains its associated data.

The optional integer argument MAXNUM limits the results to the
first MAXNUM completions. Otherwise, all completions are
returned.

The FILTER argument sets a filter function for the
completions. If supplied, it is called for each possible
completion with two arguments: the completion, and its associated
data. If the filter function returns nil, the completion is not
included in the results, and does not count towards MAXNUM."

  (let* ((node (trie--node-find tree key))
	 (trie--heap (heap-create
			`(lambda (a b) (not (,rankfun a b)))
			(when maxnum (1+ maxnum))))
	 (accumulator
	  (lambda (node seq)
	    (let ((data (trie--node-data node)))
	      (when (or (null filter) (funcall filter seq data))
		(heap-add trie--heap (cons seq data))
		(and maxnum
		     (> (heap-size trie--heap) maxnum)
		     (heap-delete-root trie--heap)))))))
    (when node (trie--mapc accumulator (trie--mapfun tree) node key))

    (let (completions)
      (while (not (heap-empty trie--heap))
	(setq completions (cons (heap-delete-root trie--heap) completions)))
      completions)))



;;; trie.el ends here
