
;;; trie.el --- trie package


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
;; Quick Overview
;; --------------
;; A trie is a data structure used to store keys that are ordered
;; sequences of elements (vectors, lists or strings in Elisp), in such a
;; way that both storage and retrieval are reasonably space- and
;; time-efficient. But, more importantly, searching for keys that match
;; various patterns can also be done efficiently. For example, returning
;; all strings with a given prefix, or searching for keys matching a
;; pattern containing wildcards, or searching for all keys within a given
;; Lewenstein distance of given string (though the latter two are not yet
;; implemented in this package - code contributions welcome!).
;;
;; You create a ternary search tree using `trie-create', create an
;; association using `trie-insert', retrieve an association using
;; `trie-lookup', find completions of a sequence using `trie-complete',
;; and map over a tree using `trie-map', `trie-mapc', `trie-mapcar', or
;; `trie-mapf'. Using `trie-stack', you can create an object that allows
;; the contents of the trie to be used like a stack; `trie-stack-pop'
;; pops elements off the stack one-by-one, whilst `trie-stack-push'
;; pushes things onto the stack.
;;
;; Note that there are two uses for a trie: as a lookup table, in which
;; case only the presence or absence of a key in the trie is significant,
;; or as an associative array, in which case each key carries some
;; associated data. Libraries for other data structure often only
;; implement lookup tables, leaving it up to you to implement an
;; associative array on top of this (by storing key+data pairs in the
;; data structure's keys, then defining a comparison function that only
;; compares the key part). However, for a trie, this would be slightly
;; less space-efficient than it needs to be, so this package does the
;; opposite: it implements associative arrays, and leaves it up to you to
;; use them as lookup tables if you so desire (with no loss of
;; space-efficiency).
;;
;;
;; Different Types of Trie
;; -----------------------
;; There are numerous ways to implement trie data structures internally,
;; each with its own trade-offs. By viewing a trie as a tree whose nodes
;; are themselves lookup tables for key elements, this package is able to
;; support all types of trie, providing there exists (or you write!) an
;; Elisp implementation of the corresponding type of lookup table. The
;; best implementation will depend on what trade-offs are appropriate for
;; your particular application. The following gives an overview of the
;; advantages and disadvantages of various types of trie. (Not all of the
;; underlying lookup tables have been implemented in Elisp yet, so using
;; some of them would require writing the missing Elisp package!)
;;
;; One of the most effective all-round implementations of a trie is a
;; ternary search tree, which can be viewed as a tree of binary trees. If
;; basic binary search trees are used for the nodes of the trie, we get a
;; basic ternary search tree. If self-balancing binary trees are used
;; (e.g. AVL or red-black trees), we get a self-balancing ternary search
;; tree. If splay trees are used, we get yet another self-organising
;; variant of a ternary search tree. All ternary search trees have, in
;; common, good space-efficiency. The time-efficiencies for the various
;; trie operations are also good, assuming the underlying binary trees
;; are balanced. Under that assumption, all variants of ternary search
;; trees described below have the same asymptotic time-complexity for all
;; trie operations.
;;
;; Self-balancing trees ensure the underlying binary trees are always
;; close to perfectly balanced, with the usual trade-offs between the
;; different the types of self-balancing binary tree: AVL trees are
;; slightly more efficient for lookup operations than red-black trees,
;; but are slightly less efficienct for insertion operations, and less
;; efficient for deletion operations. Splay trees give good average-case
;; complexity and are simpler to implement than AVL or red-black trees
;; (which can mean they're faster in practice!), at the expense of poor
;; worst-case complexity.
;;
;; If your tries are going to be static (i.e. created once and rarely
;; modified), then using perfectly balanced binary search trees might be
;; more appropriate. Perfectly balancing the binary trees is very
;; inefficient, but it only has to be when the trie is first created or
;; modified. Lookup operations will then be as efficient as possible for
;; ternary search trees, and the implementation will be much simpler (so
;; probably faster) than a self-balancing tree, without the space and
;; time overhead required to keep track of rebalancing.
;;
;; On the other hand, adding data to a binary search tree in a random
;; order usually results in a reasonably balanced tree. If this is the
;; likely scenario, using a basic binary tree without bothering to
;; balance it at all might be quite efficient, and, being even simpler to
;; implement, could be faster overall.
;;
;; A digital trie is a different implementation of a trie, which can be
;; viewed as a tree of arrays, and has different space- and
;; time-complexities than a ternary search tree. Roughly speaking, a
;; digital trie has worse space-complexity, but better
;; time-complexity. Using hash tables instead of arrays for the nodes
;; gives something similar to a digital trie, potentially with better
;; space-complexity and the same amortised time-complexity, but at the
;; expense of occasional significant inefficiency when inserting and
;; deleting (whenever the hash table has to be resized). Indeed, an array
;; can be viewed as a perfect hash table, but as such it requires the
;; number of possible values to be known in advance.
;;
;; Finally, if you really need optimal efficiency from your trie, you
;; could even write a custom type of underlying lookup table, optimised
;; for your specific needs.
;;
;;
;; This package uses the AVL tree package avl-tree.el and the heap
;; package heap.el.


;;; Change Log:
;;
;; Version 0.1
;; * Initial release (complete rewrite from scratch of tstree.el!)
;; * Ternary search trees are now implemented as a tree of avl trees, which
;;   has numerous advantages: self-balancing trees guarantee O(log n)
;;   complexity regardless of how the tree is built; deletion is now done
;;   properly.
;; * unlike tstree.el, trie.el is general enough to implement all sorts
;;   of tries, not just ternary search trees (though these remain the
;;   default).
;; * Up to "tstree"->"trie" renaming, many functions are drop-in
;;   replacements for tstree.el functions. However, insertion and rank
;;   functions are no longer stored in the data structure, so
;;   corresponidng arguments are no longer optional. A single
;;   `trie-complete' function now deals with sorting completions in both
;;   lexical or arbitrary order, the ranking function being passed as an
;;   optional argument in the latter case. And functions can no longer
;;   operate over multiple data structures at once; i.e. they no longer
;;   accept lists of trees as arguments. (These features belong in higher
;;   level libraries, and the efficiency loss is negligible.)
;; * `trie-wildcard-search' implements efficient shell-glob-like wildcard
;;   searches of tries!



;;; Code:

(eval-when-compile (require 'cl))
(require 'avl-tree)
(require 'heap)



;;; ================================================================
;;;                Setup pre-defined trie types

;; --- avl-tree ---
(put 'avl :trie-createfun (lambda (cmpfun seq) (avl-tree-create cmpfun)))
(put 'avl :trie-insertfun 'avl-tree-enter)
(put 'avl :trie-deletefun 'avl-tree-delete)
(put 'avl :trie-lookupfun 'avl-tree-member)
(put 'avl :trie-mapfun 'avl-tree-mapc)
(put 'avl :trie-emptyfun 'avl-tree-empty)
(put 'avl :trie-stack-createfun 'avl-tree-stack)
(put 'avl :trie-stack-popfun 'avl-tree-stack-pop)
(put 'avl :trie-stack-emptyfun 'avl-tree-stack-empty-p)
(put 'avl :trie-transform-for-print 'trie--avl-transform-for-print)
(put 'avl :trie-transform-from-read 'trie--avl-transform-from-read)



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


(defun trie--position (item list)
  "Find the first occurrence of ITEM in LIST.
Return the index of the matching item, or nil of not found.
Comparison is done with 'equal."
  (let (el (i 0))
    (catch 'found
      (while (setq el (nth i list))
        (when (equal item el) (throw 'found i))
        (setq i (1+ i))
        nil))))


(defsubst trie--seq-append (seq el)
  "Append EL to the end of sequence SEQ."
  (cond
   ((stringp seq) (concat seq (string el)))
   ((vectorp seq) (vconcat seq (vector el)))
   ((listp seq)	  (append seq (list el)))))


(defsubst trie--seq-concat (seq &rest sequences)
  "Concatenate SEQ and SEQUENCES, and make the result the same
type of sequence as SEQ."
  (cond
   ((stringp seq) (apply 'concat  seq sequences))
   ((vectorp seq) (apply 'vconcat seq sequences))
   ((listp seq)	  (apply 'append  seq sequences))))



;;; ================================================================
;;;     Internal functions only for use within the trie package


;;; ----------------------------------------------------------------
;;;           Functions and macros for handling a trie.

;; symbol used to denote a trie leaf node
(defconst trie--terminator '--trie--terminator)

(defstruct
  (trie-
   :named
   (:constructor nil)
   (:constructor trie--create
		 (comparison-function &optional (type 'avl)
		  &aux
		  (createfun
		   (or (get type :trie-createfun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (insertfun
		   (or (get type :trie-insertfun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (deletefun
		   (or (get type :trie-deletefun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (lookupfun
		   (or (get type :trie-lookupfun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (mapfun
		   (or (get type :trie-mapfun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (emptyfun
		   (or (get type :trie-emptyfun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (stack-createfun
		   (or (get type :trie-stack-createfun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (stack-popfun
		   (or (get type :trie-stack-popfun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (stack-emptyfun
		   (or (get type :trie-stack-emptyfun)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (transform-for-print
		   (or (get type :trie-transform-for-print)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (transform-from-read
		   (or (get type :trie-transform-from-read)
		       (error "trie--create: unknown trie TYPE, %s" type)))
		  (cmpfun (trie--wrap-cmpfun comparison-function))
		  (root (trie--node-create-root createfun cmpfun))
		  ))
   (:constructor trie--create-custom
		 (comparison-function
		  &key
		  (createfun 'avl-tree-create-bare)
		  (insertfun 'avl-tree-enter)
		  (deletefun 'avl-tree-delete)
		  (lookupfun 'avl-tree-member)
		  (mapfun 'avl-tree-mapc)
		  (emptyfun 'avl-tree-empty)
		  (stack-createfun 'avl-tree-stack)
		  (stack-popfun 'avl-tree-stack-pop)
		  (stack-emptyfun 'avl-tree-stack-empty-p)
		  (transform-for-print 'trie--avl-transform-for-print)
		  (transform-from-read 'trie--avl-transform-from-read)
		  &aux
		  (cmpfun (trie--wrap-cmpfun comparison-function))
		  (root (trie--node-create-root createfun cmpfun))
		  ))
   (:copier nil))
  root comparison-function cmpfun
  createfun insertfun deletefun lookupfun mapfun emptyfun
  stack-createfun stack-popfun stack-emptyfun
  transform-for-print transform-from-read print-form)


(defun trie--wrap-cmpfun (cmpfun)
  ;; wrap CMPFUN for use in a subtree
  `(lambda (a b)
     (setq a (trie--node-split a)
	   b (trie--node-split b))
     (cond ((eq a trie--terminator)
	    (if (eq b trie--terminator) nil t))
	   ((eq b trie--terminator) nil)
	   (t (,cmpfun a b)))))


;;; ----------------------------------------------------------------
;;;          Functions and macros for handling a trie node.

(defstruct
  (trie--node
   (:type vector)
   (:constructor nil)
   (:constructor trie--node-create
		 (split seq trie
		  &aux (subtree (funcall (trie--createfun trie)
					 (trie--cmpfun trie) seq))))
   (:constructor trie--node-create-data
		 (data &aux (split trie--terminator) (subtree data)))
   (:constructor trie--node-create-dummy
		 (split &aux (subtree nil)))
   (:constructor trie--node-create-root
		 (createfun cmpfun
		  &aux
		  (split nil)
		  (subtree (funcall createfun cmpfun []))))
   (:copier nil))
   split subtree)

;; data is stored in the subtree cell of a terminal node
(defalias 'trie--node-data 'trie--node-subtree)

(defsetf trie--node-data (node) `(setf (trie--node-subtree ,node)))

(defmacro trie--node-data-p (node)
  ;; Return t if NODE is a data node, nil otherwise.
  `(eq (trie--node-split ,node) trie--terminator))

(defmacro trie--node-p (node)
  ;; Return t if NODE is a TRIE trie--node, nil otherwise.
  ;; Have to define this ourselves, because we created a defstruct without any
  ;; identifying tags (i.e. (:type vector)) for efficiency, but this means we
  ;; can only perform a rudimentary and very unreliable test.
  `(and (vectorp ,node) (= (length ,node) 2)))


(defun trie--node-find (node seq lookupfun)
  ;; Returns the node below NODE corresponding to SEQ, or nil if none found.
  (let ((len (length seq))
	(i -1))
    ;; descend trie until we find SEQ or run out of trie
    (while (and node (< (incf i) len))
      (setq node
	    (funcall lookupfun
		     (trie--node-subtree node)
		     (trie--node-create-dummy (elt seq i))
		     nil)))
    node))


(defmacro trie--find-data-node (node lookupfun)
  ;; Return data node from NODE's subtree, or nil if NODE has no data node in
  ;; its subtree.
  `(funcall ,lookupfun
	    (trie--node-subtree ,node)
	    (trie--node-create-dummy trie--terminator)
	    nil))


(defmacro trie--find-data (node lookupfun)
  ;; Return data associated with sequence corresponding to NODE, or nil if
  ;; sequence has no associated data.
  `(let ((node (trie--find-data-node ,node ,lookupfun)))
     (when node (trie--node-data node))))



;;; ----------------------------------------------------------------
;;;              print/read transformation functions

(defun trie-transform-for-print (trie)
  "Transform TRIE to print form."
  (when (trie--transform-for-print trie)
    (if (trie--print-form trie)
	(warn "Trie has already been transformed to print-form")
      (funcall (trie--transform-for-print trie) trie)
      (setf (trie--print-form trie) t))))


(defun trie-transform-from-read (trie)
  "Transform TRIE from print form."
  (when (trie--transform-from-read trie)
    (if (not (trie--print-form trie))
	(warn "Trie is not in print-form")
      (funcall (trie--transform-from-read trie) trie)
      (setf (trie--print-form trie) nil))))


(defmacro trie-transform-from-read-warn (trie)
  "Transform TRIE from print form, with warning."
  `(when (trie--print-form ,trie)
     (warn (concat "Attempt to operate on trie in print-form; converting to\
 normal form"))
     (trie-transform-from-read ,trie)))


(defun trie--avl-transform-for-print (trie)
  ;; transform avl-tree based TRIE to print form.
  (trie-mapc-internal
   (lambda (avl seq) (setf (avl-tree--cmpfun avl) nil))
   trie))


(defun trie--avl-transform-from-read (trie)
  ;; transform avl-tree based TRIE from print form."
  (let ((--trie-avl-transform--cmpfun (trie--cmpfun trie)))
    (trie-mapc-internal
     (lambda (avl seq)
       (setf (avl-tree--cmpfun avl) --trie-avl-transform--cmpfun))
     trie)))




;;; ================================================================
;;;        The public functions which operate on tries.

(defalias 'trie-create 'trie--create
  "Return a new trie that uses comparison function COMPARISON-FUNCTION.

A trie stores sequences (strings, vectors or lists) along with
associated data. COMPARISON-FUNCTEION should accept two
arguments, each being an element of such a sequence, and return t
if the first is strictly smaller than the second.

The optional argument TYPE specifies the type of trie to
create. However, the only one that is currently implemented is
the default, so this argument is useless. (See also
`trie-create-custom'.)")



(defalias 'trie-create-custom 'trie--create-custom
  "Return a new trie that uses comparison function COMPARISON-FUNCTION.

A trie stores sequences (strings, vectors or lists) along with
associated data. COMPARISON-FUNCTION should accept two arguments,
each being an element of such a sequence, and return t if the
first is strictly smaller than the second.

The remaining keyword arguments: :CREATEFUN, :INSERTFUN, :DELETEFUN,
:LOOKUPFUN, :MAPFUN, :EMPTYFUN, :STACK-CREATEFUN, :STACK-POPFUN,
:STACK-EMPTYFUN, :TRANSFORM-FOR-PRINT and :TRANSFORM-FROM-READ
determine the type of trie that is created.

CREATEFUN is called as follows:

  (CREATEFUN COMPARISON-FUNCTION SEQ)

and should return a data structure (\"ARRAY\") that can be used
as an associative array, where two elements A and B are equal if
the following is non-nil:

  (and (COMPARISON-FUNCTION b a)
       (COMPARISON-FUNCTION b a))

The SEQ argument is a vector containing the sequence that will
correspond to the newly created array in the trie. For most types
of trie, this value is ignored. It is passed to CREATEFUN only in
order to allow the creation of \"hybrid\" trie structures, in
which different types of associative array are used in different
parts of the trie. For example, the type of associative array
could be chosen based on the depth in the trie, given by \(length
SEQ\). (Note that all the other functions described below must be
able to correctly handle *any* of the types of associate array
that might be created by CREATEFUN.)

INSERTFUN, DELETEFUN, LOOKUPFUN, MAPFUN and EMPTYFUN should
insert, delete, lookup, map over, and check-if-there-exist-any
elements in an associative array. They are called as follows:

  (INSERTFUN array element &optional updatefun)
  (DELETEFUN array element &optional predicate nilflag)
  (LOOKUPFUN array element &optional nilflag)
  (MAPFUN function array &optional reverse)
  (EMPTYFUN array)

INSERTFUN should insert ELEMENT into ARRAY and return the new
element, which will be ELEMENT itself unless UPDATEFUN is
specified. In that case, if and only if an element matching
ELEMENT already exists in the associative array, INSERTFUN should
instead pass ELEMENT and the matching element as arguments to
UPDATEFUN, replace the matching element with the return value,
and return that return value.

DELETEFUN should delete the element in the associative array that
matches ELEMENT, and return the deleted element. However, if
PREDICATE is specified and a matching element exists in ARRAY,
DELETEFUN should first pass the matching element as an argument
to PREDICATE before deleting, and should only delete the element
if PREDICATE returns non-nil. DELETEFUN should return NILFLAG if
no element was deleted (either becuase no matching element was
found, or because TESTFUN returned nil).

LOOKUPFUN should return the element from the associative array
that matches ELEMENT, or NILFLAG if no matching element exists.

MAPFUN should map FUNCTION over all elements in the order defined by
COMPARISON-FUNCTION, or in reverse order if REVERSE is non-nil.


STACK-CREATEFUN, STACK-POPFUN and STACK-EMPTYFUN should allow the
associative array to be used as a stack. STACK-CREATEFUN is
called as follows:

  (STACK-CREATEFUN array)

and should return a data structure (\"STACK\") that behaves like
a sorted stack of all elements in the associative array. I.e.
successive calls to

  (STACK-POPFUN stack)

should return elements from the associative array in the order
defined by COMPARISON-FUNCTION, and

  (STACK-EMPTYFUN stack)

should return non-nil if the stack is empty, nil otherwise.

The stack functions are optional, in that all trie operations
other than the stack-related ones will work correctly. However,
any code that makes use of trie-stacks will complain if supplied
with this type of trie.


The :TRANSFORM-FOR-PRINT and :TRANSFORM-FROM-READ arguments are
optional. If supplied, they can be used to transform the trie
into a format suitable for passing to Elisp's `print'
functions (typically used to persistently store the trie by
writing it to file), and transform from that format back to the
original usable form.


Warning: to avoid nasty dynamic scoping bugs, the supplied
functions must *never* bind any variables with names commencing \"--\".")



(defalias 'trie-comparison-function 'trie--comparison-function
  "Return the comparison function for TRIE.")


(defalias 'trie-p 'trie--p
  "Return t if argument is a trie, nil otherwise.")


(defun trie-empty (trie)
  "Return t if the TRIE is empty, nil otherwise."
  (trie-transform-from-read-warn trie)
  (funcall (trie--emptyfun trie) (trie--node-subtree (trie--root trie))))


(defun trie-construct-sortfun (cmpfun &optional reverse)
  "Construct function to compare key sequences, based on a CMPFUN
that compares individual elements of the sequence. Order is
reversed if REVERSE is non-nil."
  (if reverse
      `(lambda (a b)
	 (let (cmp)
	   (catch 'compared
	     (dotimes (i (min (length a) (length b)))
	       (cond ((,cmpfun (elt b i) (elt a i)) (throw 'compared t))
		     ((,cmpfun (elt a i) (elt b i)) (throw 'compared nil))))
	     (< (length a) (length b)))))
    `(lambda (a b)
       (let (cmp)
	 (catch 'compared
	   (dotimes (i (min (length a) (length b)))
	     (cond ((,cmpfun (elt a i) (elt b i)) (throw 'compared t))
		   ((,cmpfun (elt b i) (elt a i)) (throw 'compared nil))))
	   (< (length a) (length b)))))))



;; ----------------------------------------------------------------
;;                        Inserting data

(defun trie-insert (trie key &optional data updatefun)
  "Associate DATA with KEY in TRIE.

If KEY already exists in TRIE, then DATA replaces the existing
association, unless UPDATEFUN is supplied. Note that if DATA is
*not* supplied, this means that the existing association of KEY
will be replaced by nil.

If UPDATEFUN is supplied and KEY already exists in TRIE,
UPDATEFUN is called with two arguments: DATA and the existing
association of KEY. Its return value becomes the new association
for KEY.

Returns the new association of KEY.

Note: to avoid nasty dynamic scoping bugs, UPDATEFUN must *not*
bind any variables with names commencing \"--\"."

  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)

  ;; absurd variable names are an attempt to avoid dynamic scoping bugs
  (let ((--trie-insert--updatefun updatefun)
	--trie-insert--old-node-flag
	(node (trie--root trie))
	(len (length key))
	(i -1))
    ;; Descend trie, adding nodes for non-existent elements of KEY. The
    ;; update function passed to `trie--insertfun' ensures that existing
    ;; nodes are left intact.
    (while (< (incf i) len)
      (setq --trie-insert--old-node-flag nil)
      (setq node (funcall (trie--insertfun trie)
			  (trie--node-subtree node)
			  (trie--node-create (elt key i) key trie)
			  (lambda (a b)
			    (setq --trie-insert--old-node-flag t) b))))
    ;; Create or update data node.
    (setq node (funcall (trie--insertfun trie)
			(trie--node-subtree node)
			(trie--node-create-data data)
			;; if using existing data node, wrap UPDATEFUN if
			;; any was supplied
			(when (and --trie-insert--old-node-flag
				   --trie-insert--updatefun)
			  (lambda (new old)
			    (setf (trie--node-data old)
				  (funcall --trie-insert--updatefun
					   (trie--node-data new)
					   (trie--node-data old)))
			    old))))
    (trie--node-data node)))  ; return new data



;; ----------------------------------------------------------------
;;                        Deleting data

(defun trie-delete (trie key &optional test)
  "Delete KEY and its associated data from TRIE.

If KEY was deleted, a cons cell containing KEY and its
association is returned. Returns nil if KEY does not exist in
TRIE.

If TEST is supplied, it should be a function that accepts two
arguments: the key being deleted, and its associated data. The
key will then only be deleted if TEST returns non-nil.

Note: to avoid nasty dynamic scoping bugs, TEST must *not* bind
any variables with names commencing \"--\"."
  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; set up deletion (real work is done by `trie--do-delete'
  (let (--trie-deleted--node
	(--trie-delete--key key))
    (declare (special --trie-deleted--node)
	     (special --trie-delete--key))
    (trie--do-delete (trie--root trie) key test
		     (trie--deletefun trie)
		     (trie--emptyfun trie)
		     (trie--cmpfun trie))
    (when --trie-deleted--node
      (cons key (trie--node-data --trie-deleted--node)))))


(defun trie--do-delete (node --trie--do-delete--seq
			     --trie--do-delete--test
			     --trie--do-delete--deletefun
			     --trie--do-delete--emptyfun
			     --trie--do-delete--cmpfun)
  ;; Delete --TRIE--DO-DELETE--SEQ starting from trie node NODE, and return
  ;; non-nil if we deleted a node. If --TRIE--DO-DELETE--TEST is supplied, it
  ;; is called with two arguments, the key being deleted and the associated
  ;; data, and the deletion is only carried out if it returns non-nil.

  ;; The absurd argument names are to lessen the likelihood of dynamical
  ;; scoping bugs caused by a supplied function binding a variable with the
  ;; same name as one of the arguments, which would cause a nasty bug when the
  ;; lambda's (below) are called.
  (declare (special --trie-deleted--node)
	   (special --trie-delete--key))
  ;; if --TRIE--DO-DELETE--SEQ is empty, try to delete data node and return
  ;; non-nil if we did (return value of --TRIE--DO-DELETE--DELETEFUN is the
  ;; deleted data, which is always non-nil for a trie)
  (if (= (length --trie--do-delete--seq) 0)
      (setq --trie-deleted--node
	    (funcall --trie--do-delete--deletefun
		     (trie--node-subtree node)
		     (trie--node-create-dummy trie--terminator)
		     (when --trie--do-delete--test
		       (lambda (n)
			 (funcall --trie--do-delete--test
				  --trie-delete--key (trie--node-data n))))
		     nil))
    ;; otherwise, delete on down (return value of --TRIE--DO-DELETE--DELETEFUN
    ;; is the deleted data, which is always non-nil for a trie)
    (funcall --trie--do-delete--deletefun
	     (trie--node-subtree node)
	     (trie--node-create-dummy (elt --trie--do-delete--seq 0))
	     (lambda (n)
	       (and (trie--do-delete
		     n (trie--subseq --trie--do-delete--seq 1)
		     --trie--do-delete--test
		     --trie--do-delete--deletefun
		     --trie--do-delete--emptyfun
		     --trie--do-delete--cmpfun)
		    (funcall --trie--do-delete--emptyfun
			     (trie--node-subtree n))))
	     nil)))



;; ----------------------------------------------------------------
;;                       Retrieving data

(defun trie-lookup (trie key &optional nilflag)
  "Return the data associated with KEY in the TRIE,
or nil if KEY does not exist in TRIE.

Optional argument NILFLAG specifies a value to return instead of
nil if KEY does not exist in TRIE. This allows a non-existent KEY
to be distinguished from an element with a null association. (See
also `trie-member-p', which does this for you.)"
  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; find node corresponding to key, then find data node, then return data
  (let (node)
    (or (and (setq node (trie--node-find (trie--root trie) key
					 (trie--lookupfun trie)))
	     (trie--find-data node (trie--lookupfun trie)))
	nilflag)))

(defalias 'trie-member 'trie-lookup)


(defun trie-member-p (trie key)
  "Return t if KEY exists in TRIE, nil otherwise."
  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)
  (let ((flag '(nil)))
    (not (eq flag (trie-member trie key flag)))))



;;; ----------------------------------------------------------------
;;;                      Mapping over tries

(defun trie--mapc (--trie--mapc--function --trie--mapc--mapfun
		   --trie--mapc--root --trie--mapc--seq
		   &optional --trie--mapc--reverse)
  ;; Apply TRIE--MAPC--FUNCTION to all elements in a trie beneath
  ;; TRIE--MAPC--ROOT, which should correspond to the sequence
  ;; TRIE--MAPC--SEQ. TRIE--MAPC--FUNCTION is passed two arguments: the trie
  ;; node itself and the sequence it corresponds to. It is applied in
  ;; ascending order, or descending order if TRIE--MAPC--REVERSE is non-nil.

  ;; The absurd argument names are to lessen the likelihood of dynamical
  ;; scoping bugs caused by a supplied function binding a variable with the
  ;; same name as one of the arguments.
  (funcall
   --trie--mapc--mapfun
   (lambda (--trie--mapc--node)
     ;; data node: apply function
     (if (trie--node-data-p --trie--mapc--node)
	 (funcall --trie--mapc--function --trie--mapc--node --trie--mapc--seq)
       ;; internal node: append split value to seq and keep descending
       (trie--mapc --trie--mapc--function --trie--mapc--mapfun
		   --trie--mapc--node
		   (trie--seq-append (copy-sequence --trie--mapc--seq)
				     (trie--node-split --trie--mapc--node))
		   --trie--mapc--reverse)))
   ;; --TRIE--MAPC--MAPFUN target
   (trie--node-subtree --trie--mapc--root)
   --trie--mapc--reverse))


(defun trie-mapc-internal (function trie &optional type)
  "Apply FUNCTION to all internal associative arrays within TRIE.
FUNCTION is passed two arguments: an associative array, and the
sequence it corresponds to.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to function. Defaults to
vector."
  (trie--mapc-internal function (trie--mapfun trie) (trie--root trie)
		       (cond ((eq type 'string) "")
			     ((eq type 'lisp) ())
			     (t []))))


(defun trie--mapc-internal (--trie--mapc-internal--function
			     --trie--mapc-internal--mapfun
			     --trie--mapc-internal--root
			     --trie--mapc-internal--seq)
  (funcall
   --trie--mapc-internal--mapfun
   (lambda (--trie--mapc-internal--node)
     ;; data node
     (unless (trie--node-data-p --trie--mapc-internal--node)
       (funcall --trie--mapc-internal--function
		(trie--node-subtree --trie--mapc-internal--node)
		--trie--mapc-internal--seq)
       (trie--mapc-internal
	--trie--mapc-internal--function
	--trie--mapc-internal--mapfun
	--trie--mapc-internal--node
	(trie--seq-append (copy-sequence --trie--mapc-internal--seq)
			  (trie--node-split --trie--mapc-internal--node)))))
   (trie--node-subtree --trie--mapc-internal--root)))


(defun trie-map (function trie &optional type reverse)
  "Modify all elements in TRIE by applying FUNCTION to them.

FUNCTION should take two arguments: a sequence stored in the trie
and its associated data. Its return value replaces the existing
data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

FUNCTION is applied in ascending order, or descending order if
REVERSE is non-nil.

Note: to avoid nasty dynamic scoping bugs, FUNCTION must *not*
bind any variables with names commencing \"--\"."
  ;; convert from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; map FUNCTION over TRIE
  (let ((--trie-map--function function)) ; try to avoid dynamic scoping bugs
    (trie--mapc
     (lambda (node seq)
       (setf (trie--node-data node)
	     (funcall --trie-map--function seq (trie--node-data node))))
     (trie--mapfun trie)
     (trie--root trie)
     (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
     reverse)))


(defun trie-mapc (function trie &optional type reverse)
  "Apply FUNCTION to all elements in TRIE for side effect only.

FUNCTION should take two arguments: a sequence stored in the trie
and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

FUNCTION is applied in ascending order, or descending order if
REVERSE is non-nil.

Note: to avoid nasty dynamic scoping bugs, FUNCTION must *not*
bind any variables with names commencing \"--\"."
  ;; convert from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; map FUNCTION over TRIE
  (let ((--trie-mapc--function function)) ; try to avoid dynamic scoping bugs
    (trie--mapc
     (lambda (node seq)
       (funcall --trie-mapc--function seq (trie--node-data node)))
     (trie--mapfun trie)
     (trie--root trie)
     (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
     reverse)))


(defun trie-mapf (function combinator trie &optional type reverse)
  "Apply FUNCTION to all elements in TRIE, and combine the results
using COMBINATOR.

FUNCTION should take two arguments: a sequence stored in the
trie, and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string; defaults to vector) sets the type of sequence passed to
FUNCTION. If TYPE is 'string, it must be possible to apply the
function `string' to the individual elements of key sequences
stored in TRIE.

The FUNCTION is applied and the results combined in ascending
order, or descending order if REVERSE is non-nil.

Note: to avoid nasty dynamic scoping bugs, FUNCTION and
COMBINATOR must *not* bind any variables with names
commencing \"--\"."
  ;; convert from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; map FUNCTION over TRIE, combining results with COMBINATOR
  (let ((--trie-mapf--function function) ; try to avoid dynamic scoping bugs
	--trie-mapf--accumulate)
    (trie--mapc
     (lambda (node seq)
       (setq --trie-mapf--accumulate
	     (funcall combinator
		      (funcall --trie-mapf--function
			       seq (trie--node-data node))
		      --trie-mapf--accumulate)))
     (trie--mapfun trie)
     (trie--root trie)
     (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
     reverse)
    --trie-mapf--accumulate))


(defun trie-mapcar (function trie &optional type reverse)
  "Apply FUNCTION to all elements in TRIE, and make a list of the results.

FUNCTION should take two arguments: a sequence stored in the trie
and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

The FUNCTION is applied and the list constructed in ascending
order, or descending order if REVERSE is non-nil.

Note that if you don't care about the order in which FUNCTION is
applied, just that the resulting list is in the correct order,
then

  (trie-mapf function 'cons trie type (not reverse))

is more efficient.

Note: to avoid nasty dynamic scoping bugs, FUNCTION must *not*
bind any variables with names commencing \"--\"."
  ;; convert from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; map FUNCTION over TRIE and accumulate in a list
  (nreverse (trie-mapf function 'cons trie type reverse)))



;;; ----------------------------------------------------------------
;;;                    Using tries as stacks

(defstruct (trie--stack
	    (:constructor nil)
	    (:constructor
	     trie--stack-create
	     (trie
	      &optional
	      (type 'vector)
	      reverse
	      &aux
	      (comparison-function (trie--comparison-function trie))
	      (lookupfun (trie--lookupfun trie))
	      (stack-createfun (trie--stack-createfun trie))
	      (stack-popfun (trie--stack-popfun trie))
	      (stack-emptyfun (trie--stack-emptyfun trie))
	      (repopulatefun 'trie--stack-repopulate)
	      (store
	       (if (trie-empty trie)
		   nil
		 (trie--stack-repopulate
		  (list (cons
			 (cond ((eq type 'list) ())
			       ((eq type 'string) "")
			       (t []))
			 (funcall stack-createfun
				  (trie--node-subtree (trie--root trie))
				  reverse)))
		  reverse
		  comparison-function lookupfun
		  stack-createfun stack-popfun stack-emptyfun)))
	      (pushed '())
	      ))
	    (:constructor
	     trie--completion-stack-create
	     (trie prefix
	      &optional
	      reverse
	      &aux
	      (comparison-function (trie--comparison-function trie))
	      (lookupfun (trie--lookupfun trie))
	      (stack-createfun (trie--stack-createfun trie))
	      (stack-popfun (trie--stack-popfun trie))
	      (stack-emptyfun (trie--stack-emptyfun trie))
	      (repopulatefun 'trie--stack-repopulate)
	      (store (trie--completion-stack-construct-store
		      trie prefix reverse))
	      (pushed '())
	      ))
	    (:constructor
	     trie--wildcard-stack-create
	     (trie pattern
	      &optional
	      reverse
	      &aux
	      (comparison-function (trie--comparison-function trie))
	      (lookupfun (trie--lookupfun trie))
	      (stack-createfun (trie--stack-createfun trie))
	      (stack-popfun (trie--stack-popfun trie))
	      (stack-emptyfun (trie--stack-emptyfun trie))
	      (repopulatefun 'trie--wildcard-stack-repopulate)
	      (store (trie--wildcard-stack-construct-store
		      trie pattern reverse))
	      (pushed '())
	      ))
	    (:copier nil))
  reverse comparison-function lookupfun
  stack-createfun stack-popfun stack-emptyfun
  repopulatefun store pushed)


(defun trie-stack (trie &optional type reverse)
  "Return an object that allows TRIE to be accessed as if it were a stack.

The stack is sorted in \"lexical\" order, i.e. the order defined
by the trie's comparison function, or in reverse order if REVERSE
is non-nil. Calling `trie-stack-pop' pops the top element (a key
and its associated data) from the stack.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence used for the keys.

Note that any modification to TRIE *immediately* invalidates all
trie-stacks created before the modification (in particular,
calling `trie-stack-pop' will give unpredictable results).

Operations on trie-stacks are significantly more efficient than
constructing a real stack from the trie and using standard stack
functions. As such, they can be useful in implementing efficient
algorithms on tries. However, in cases where mapping functions
`trie-mapc', `trie-mapcar' or `trie-mapf' would be sufficient, it
is better to use one of those instead."
  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; if stack functions aren't defined for trie type, throw error
  (if (not (functionp (trie--stack-createfun trie)))
      (error "Trie type does not support stack operations")
    ;; otherwise, create and initialise a stack
    (trie--stack-create trie type reverse)))


(defun trie-stack-pop (trie-stack)
  "Pop the first element from TRIE-STACK.
Returns nil if the stack is empty."
  ;; if elements have been pushed onto the stack, pop those first
  (if (trie--stack-pushed trie-stack)
      (pop (trie--stack-pushed trie-stack))
    ;; otherwise, pop first element from trie-stack and repopulate it
    (let ((first (pop (trie--stack-store trie-stack))))
      (when first
	(setf (trie--stack-store trie-stack)
	      (funcall (trie--stack-repopulatefun trie-stack)
		       (trie--stack-store trie-stack)
		       (trie--stack-reverse trie-stack)
		       (trie--stack-comparison-function trie-stack)
		       (trie--stack-lookupfun trie-stack)
		       (trie--stack-stack-createfun trie-stack)
		       (trie--stack-stack-popfun trie-stack)
		       (trie--stack-stack-emptyfun trie-stack)))
	first))))


(defun trie-stack-push (element trie-stack)
  "Push ELEMENT onto TRIE-STACK."
  (push element (trie--stack-pushed trie-stack)))


(defun trie-stack-first (trie-stack)
  "Return the first element from TRIE-STACK, without removing it
from the stack. Returns nil if the stack is empty."
  ;; if elements have been pushed onto the stack, return first of those
  (if (trie--stack-pushed trie-stack)
      (car (trie--stack-pushed trie-stack))
    ;; otherwise, return first element from trie-stack
    (car (trie--stack-store trie-stack))))


(defalias 'trie-stack-p 'trie--stack-p
  "Return t if argument is a trie-stack, nil otherwise.")


(defun trie-stack-empty-p (trie-stack)
  "Return t if TRIE-STACK is empty, nil otherwise."
  (and (null (trie--stack-store trie-stack))
       (null (trie--stack-pushed trie-stack))))


(defun trie--stack-repopulate (store reverse
			       comparison-function lookupfun
			       stack-createfun stack-popfun stack-emptyfun)
  ;; Recursively push children of the node at the head of STORE onto the front
  ;; of STORE, until a data node is reached.

  ;; nothing to do if stack is empty
  (when store
    (let ((node (funcall stack-popfun (cdar store)))
	  (seq (caar store)))
      (when (funcall stack-emptyfun (cdar store))
	;; (pop store) here produces irritating compiler warnings
	(setq store (cdr store)))

      (while (not (trie--node-data-p node))
	(push
	 (cons (trie--seq-append seq (trie--node-split node))
	       (funcall stack-createfun (trie--node-subtree node) reverse))
	 store)
	(setq node (funcall stack-popfun (cdar store))
	      seq (caar store))
	(when (funcall stack-emptyfun (cdar store))
	  ;; (pop store) here produces irritating compiler warnings
	  (setq store (cdr store))))

      (push (cons seq (trie--node-data node)) store))))



;; ----------------------------------------------------------------
;;                   Advanced query-building macros

;; Implementation Note
;; -------------------
;; For queries ranked in anything other than lexical order, we use a partial
;; heap-sort to find the k=MAXNUM highest ranked matches among the n possibile
;; matches. This has worst-case time complexity O(n log k), and is both simple
;; and elegant. An optimal algorithm (e.g. partial quick-sort where the
;; irrelevant partition is discarded at each step) would have complexity O(n +
;; k log k), but is probably not worth the extra coding effort, and would have
;; worse space complexity unless coded to work "in-place", which would be
;; highly non-trivial. (I haven't done any benchmarking, though, so feel free
;; to do so and let me know the results!)

(defmacro trie--construct-accumulator (maxnum filter)
  ;; Does what it says on the tin! | sed -e 's/on/in/' -e 's/tin/macro name/'
  `(cond
    ((and ,filter ,maxnum)
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (when (funcall ,filter seq data)
	   (aset trie--accumulate 0
		 (cons (cons seq data)
		       (aref trie--accumulate 0)))
	   (and (>= (length (aref trie--accumulate 0)) ,maxnum)
		(throw 'trie-complete--done nil))))))
    ((and (not ,filter) ,maxnum)
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (aset trie--accumulate 0
	       (cons (cons seq data)
		     (aref trie--accumulate 0)))
	 (and (>= (length (aref trie--accumulate 0)) ,maxnum)
	      (throw 'trie-complete--done nil)))))
    ((and ,filter (not ,maxnum))
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (when (funcall ,filter seq data)
	   (aset trie--accumulate 0
		 (cons (cons seq data)
		       (aref trie--accumulate 0)))))))
    ((and (not ,filter) (not ,maxnum))
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (aset trie--accumulate 0
	       (cons (cons seq data)
		     (aref trie--accumulate 0))))))))


(defmacro trie--construct-ranked-accumulator (maxnum filter)
  ;; Does what it says on the tin! | sed -e 's/on/in/' -e 's/tin/macro name/'
  `(cond
    ((and ,filter ,maxnum)
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (when (funcall ,filter seq data)
	   (heap-add trie--accumulate (cons seq data))
	   (and (> (heap-size trie--accumulate) ,maxnum)
		(heap-delete-root trie--accumulate))))))
    ((and ,filter (not ,maxnum))
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (when (funcall ,filter seq data)
	   (heap-add trie--accumulate (cons seq data))))))
    ((and (not ,filter) ,maxnum)
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (heap-add trie--accumulate (cons seq data))
	 (and (> (heap-size trie--accumulate) ,maxnum)
	      (heap-delete-root trie--accumulate)))))
    ((and (not ,filter) (not ,maxnum))
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (heap-add trie--accumulate (cons seq data)))))))



(defmacro trie--accumulate-results
  (rankfun maxnum reverse filter accfun duplicates &rest body)
  ;; Accumulate results of running BODY code, and return them in appropriate
  ;; order. BODY should call ACCFUN to accumulate a result, passing it two
  ;; arguments: a trie data node, and the corresponding sequence. A non-null
  ;; DUPLICATES flag signals that the accumulated results might contain
  ;; duplicates, which should be deleted. Note that DUPLICATES is ignored if
  ;; RANKFUN is null. The other arguments should be passed straight through
  ;; from the query function.

  ;; rename RANKFUN to help avoid dynamic-scoping bugs
  `(let* ((--trie-accumulate--rankfun ,rankfun)
	  ;; construct structure in which to accumulate results
	  (trie--accumulate
	   (if ,rankfun
	       (heap-create  ; heap order is inverse of rank order
		(if ,reverse
		    (lambda (a b)
		      (funcall --trie-accumulate--rankfun a b))
		  (lambda (a b)
		    (not (funcall --trie-accumulate--rankfun a b))))
		(when ,maxnum (1+ ,maxnum)))
	     (make-vector 1 nil)))
	  ;; construct function to accumulate completions
	  (,accfun
	   (if ,rankfun
	       (trie--construct-ranked-accumulator ,maxnum ,filter)
	     (trie--construct-accumulator ,maxnum ,filter))))

     ;; accumulate results
     (catch 'trie-complete--done ,@body)

     ;; return list of completions
     (cond
      ;; extract completions from heap for ranked query
      (,rankfun
       (let (completions)
	 ;; check for and delete duplicates if flag is set
	 (if ,duplicates
	     (while (not (heap-empty trie--accumulate))
	       (if (equal (car (heap-root trie--accumulate))
			  (caar completions))
		   (heap-delete-root trie--accumulate)
		 (push (heap-delete-root trie--accumulate) completions)))
	   ;; skip duplicate checking if flag is not set
	   (while (not (heap-empty trie--accumulate))
	     (push (heap-delete-root trie--accumulate) completions)))
	 completions))

      ;; for lexical query, reverse result list if MAXNUM supplied
      (,maxnum (nreverse (aref trie--accumulate 0)))
      ;; otherwise, just return list
      (t (aref trie--accumulate 0)))))



;; ----------------------------------------------------------------
;;                          Completing

(defun trie-complete (trie prefix &optional rankfun maxnum reverse filter)
  "Return an alist containing all completions of PREFIX in TRIE
along with their associated data, in the order defined by
RANKFUN, defaulting to \"lexical\" order (i.e. the order defined
by the trie's comparison function). If REVERSE is non-nil, the
completions are sorted in the reverse order. If no completions
are found, return nil.

PREFIX must be a sequence (vector, list or string) containing
elements of the type used to reference data in the trie. (If
PREFIX is a string, it must be possible to apply `string' to
individual elements of the sequences stored in the trie.) The
completions returned in the alist will be sequences of the same
type as KEY. If PREFIX is a list of sequences, completions of all
sequences in the list are included in the returned alist. All
sequences in the list must be of the same type.

The optional integer argument MAXNUM limits the results to the
first MAXNUM completions. Otherwise, all completions are
returned.

If specified, RANKFUN must accept two arguments, both cons
cells. The car contains a sequence from the trie (of the same
type as PREFIX), the cdr contains its associated data. It should
return non-nil if first argument is ranked strictly higher than
the second, nil otherwise.

The FILTER argument sets a filter function for the
completions. If supplied, it is called for each possible
completion with two arguments: the completion, and its associated
data. If the filter function returns nil, the completion is not
included in the results, and does not count towards MAXNUM."

  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; wrap prefix in a list if necessary
  ;; FIXME: the test for a list of prefixes, below, will fail if the PREFIX
  ;;        sequence is a list, and the elements of PREFIX are themselves
  ;;        lists (there might be no easy way to fully fix this...)
  (if (or (atom prefix)
	  (and (listp prefix) (not (sequencep (car prefix)))))
      (setq prefix (list prefix))
    ;; sort list of prefixes if sorting completions lexically
    (when (null rankfun)
      (setq prefix
	    (sort prefix (trie-construct-sortfun
			  (trie--comparison-function trie))))))

  ;; accumulate completions
  (let (node)
    (trie--accumulate-results
     rankfun maxnum reverse filter accumulator nil
     (mapc (lambda (pfx)
	     (setq node (trie--node-find (trie--root trie) pfx
					 (trie--lookupfun trie)))
	     (when node
	       (trie--mapc accumulator (trie--mapfun trie) node pfx
			   (if maxnum reverse (not reverse)))))
	   prefix))
    ))


(defun trie-complete-stack (trie prefix &optional reverse)
  "Return an object that allows completions of PREFIX to be accessed
as if they were a stack.

The stack is sorted in \"lexical\" order, i.e. the order defined
by TRIE's comparison function, or in reverse order if REVERSE is
non-nil. Calling `trie-stack-pop' pops the top element (a key and
its associated data) from the stack.

PREFIX must be a sequence (vector, list or string) that forms the
initial part of a TRIE key. (If PREFIX is a string, it must be
possible to apply `string' to individual elements of TRIE keys.)
The completions returned in the alist will be sequences of the
same type as KEY. If PREFIX is a list of sequences, completions
of all sequences in the list are included in the stack. All
sequences in the list must be of the same type.

Note that any modification to TRIE *immediately* invalidates all
trie-stacks created before the modification (in particular,
calling `trie-stack-pop' will give unpredictable results).

Operations on trie-stacks are significantly more efficient than
constructing a real stack from completions of PREFIX in TRIE and
using standard stack functions. As such, they can be useful in
implementing efficient algorithms on tries. However, in cases
where `trie-complete' or `trie-complete-ordered' is sufficient,
it is better to use one of those instead."
  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; if stack functions aren't defined for trie type, throw error
  (if (not (functionp (trie--stack-createfun trie)))
      (error "Trie type does not support stack operations")
    ;; otherwise, create and initialise a stack
    (trie--completion-stack-create trie prefix reverse)))


(defun trie--completion-stack-construct-store (trie prefix reverse)
  ;; Construct store for completion stack based on TRIE.
  (let (store node)
    (if (or (atom prefix)
	    (and (listp prefix)
		 (not (sequencep (car prefix)))))
	(setq prefix (list prefix))
      (setq prefix
	    (sort prefix
		  (trie-construct-sortfun
		   (trie--comparison-function trie)
		   (not reverse)))))
    (dolist (pfx prefix)
      (when (setq node (trie--node-find (trie--root trie) pfx
					(trie--lookupfun trie)))
	(push (cons pfx (funcall (trie--stack-createfun trie)
				 (trie--node-subtree node)
				 reverse))
	      store)))
    (trie--stack-repopulate
     store reverse
     (trie--comparison-function trie)
     (trie--lookupfun trie)
     (trie--stack-createfun trie)
     (trie--stack-popfun trie)
     (trie--stack-emptyfun trie))))



;; ----------------------------------------------------------------
;;                        Wildcard search

(defmacro trie--wildcard-literal-p (el) `(vectorp ,el))

(defmacro trie--wildcard-*-p (el) `(eq ,el ?*))

(defmacro trie--wildcard-?-p (el) `(eq ,el ??))

(defmacro trie--wildcard-char-alt-p (el)
  `(and (listp ,el)
	(or (not (eq (car (last ,el)) ?^))
	    (= (length ,el) 1))))

(defmacro trie--wildcard-neg-char-alt-p (el)
  `(and (listp ,el)
	(eq (car (last ,el)) ?^)
	(not (= (length ,el) 1))))



(defun trie-wildcard-search (trie pattern
				  &optional rankfun maxnum reverse filter)
  "Return an alist containing all matches for PATTERN in TRIE
along with their associated data, in the order defined by
RANKFUN, defaulting to \"lexical\" order (i.e. the order defined
by the trie's comparison function). If REVERSE is non-nil, the
completions are sorted in the reverse order. If no completions
are found, return nil.

PATTERN must be a sequence (vector, list or string) containing
either elements of the type used to reference data in the trie,
or any the characters `*', `?', `[', `]', `^' or `\\'. The
meaning and syntax of these special characters follows shell-glob
syntax:

  *  wildcard
    Matches zero or more characters.

  ?  wildcard
    Matches a single character.

  [...]  character alternative
    Matches any of the listed characters.

  [^...]  negated character alternative
    Matches any character *other* then those listed.

  []..]  character alternative including `]'
    Matches any of the listed characters, including `]'.

  \\  quote literal
    Causes the next element of the pattern sequence to be treated
    literally; special characters lose their special meaning, for
    anything else it has no effect.

To include a `]' in a character alternative, place it immediately
after the opening `['. To include a literal `\\', quote it with
another `\\' (remember that `\\' also has to be quoted within
elisp strings, so as a string this would be \"\\\\\\\\\"). The
above syntax descriptions are written in terms of strings, but
the special characters can be used in *any* sequence
type. E.g. the character alternative \"[abc]\" would be \(?[ ?a
?b ?c ?]\) as a list, or [?[ ?a ?b ?c ?]] as a vector. The
\"characters\" in the alternative can of course be any data type
that might be stored in the trie, not just actual characters.

If PATTERN is a string, it must be possible to apply `string' to
individual elements of the sequences stored in the trie. The
matches returned in the alist will be sequences of the same type
as KEY. If PATTERN is a list of pattern sequences, matches for
all patterns in the list are included in the returned alist. All
sequences in the list must be of the same type.

The optional integer argument MAXNUM limits the results to the
first MAXNUM matches. Otherwise, all matches are returned.

If specified, RANKFUN must accept two arguments, both cons
cells. The car contains a sequence from the trie (of the same
type as PREFIX), the cdr contains its associated data. It should
return non-nil if first argument is ranked strictly higher than
the second, nil otherwise.

The FILTER argument sets a filter function for the matches. If
supplied, it is called for each possible match with two
arguments: the matching key, and its associated data. If the
filter function returns nil, the match is not included in the
results, and does not count towards MAXNUM.


Efficiency concerns:

Wildcard searches on tries are very efficient compared to similar
searches on other data structures. However, some wildcard
patterns are inherently time-consuming to match, especially those
containing `*' wildcards. As a general rule, patterns containing
a `*' wildcard will be slower the closer the `*' is to the
beginning of the pattern, and patterns containing multiple `*'
wildcards can be very slow indeed."

  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; wrap prefix in a list if necessary
  ;; FIXME: the test for a list of prefixes, below, will fail if the PREFIX
  ;;        sequence is a list, and the elements of PREFIX are themselves
  ;;        lists (there might be no easy way to fully fix this...)
  (if (or (atom pattern)
	  (and (listp pattern) (not (sequencep (car pattern)))))
      (setq pattern (list pattern))
    ;; sort list of patterns if sorting completions lexically
    (when (null rankfun)
      (setq pattern
	    (sort pattern (trie-construct-sortfun
			  (trie--comparison-function trie))))))

  ;; accumulate pattern matches
  (declare (special accumulator))
  (let (duplicates
	(initseq (cond ((stringp (car pattern)) "")
		       ((listp (car pattern)) ())
		       (t []))))
    ;; check for * wildcards in pattern
    (setq pattern
	  (mapcar (lambda (pat)
		    ;; convert pattern to list
		    (setq pat (append pat nil))
		    (let ((pos (trie--position ?* pat)))
		      ;; if *'s appear in middle, have to sort manually
		      (when (and (null rankfun) pos
				 (not (= pos (1- (length pat)))))
			(setq rankfun
			      `(lambda (a b)
				 (,(trie-construct-sortfun
				    (trie--comparison-function trie)
				    reverse)
				  (car a) (car b)))))
		      ;; if pattern contains multiple *'s, might get dups
		      (when (and pos
				 (trie--position
				  ?* (trie--subseq pat (1+ pos))))
			(setq duplicates t)))
		    ;; return pattern converted to list
		    pat)
		  pattern))

    (trie--accumulate-results
     rankfun maxnum reverse filter accumulator duplicates
     (mapc (lambda (pat)
	     (trie--do-wildcard-search
	      (trie--root trie)
	      initseq pat rankfun maxnum reverse
	      (trie--comparison-function trie)
	      (trie--lookupfun trie)
	      (trie--mapfun trie)))
	   pattern))))



(defun trie--do-wildcard-search
  (node seq pattern rankfun maxnum reverse
	comparison-function lookupfun mapfun)
  ;; Perform wildcard search for PATTERN starting at NODE which corresponds to
  ;; SEQ. RANKFUN, MAXNUM and REVERSE should be passed through from query
  ;; function, COMPARISON-FUNCTION, LOOKUPFUN and MAPFUN should be
  ;; corresponding trie functions (note that COMPARISON-FUNCTION should be the
  ;; trie--comparison-function, *not* the trie--cmpfun)
  (declare (special accumulator))

  ;; if pattern is null, accumulate data from current node
  (if (null pattern)
      (when (setq node (trie--find-data-node node lookupfun))
	(funcall accumulator node seq))

    ;; otherwise, extract first pattern element and act on it
    (setq pattern (trie--wildcard-parse-pattern pattern))
    (let ((el (car pattern)))
      (setq pattern (cdr pattern))

      (cond
       ;; literal string: descend to corresponding node
       ((trie--wildcard-literal-p el)
	(when (setq node (trie--node-find node el lookupfun))
	  (trie--do-wildcard-search node (trie--seq-concat seq el)
				    pattern rankfun maxnum reverse
				    comparison-function lookupfun mapfun)))

       ;; terminal *: accumulate everything below current node
       ((and (null pattern) (trie--wildcard-*-p el))
	(trie--mapc accumulator mapfun node seq
		    (if maxnum reverse (not reverse))))

       ;; * wildcard: map over all nodes immediately below current one, with
       ;;             and without using up the *
       ((trie--wildcard-*-p el)
	(funcall mapfun
		 (lambda (node)
		   ;; skip data nodes (terminal * dealt with above)
		   (unless (trie--node-data-p node)
		     ;; using up *
		     (trie--do-wildcard-search
		      node (trie--seq-append seq (trie--node-split node))
		      pattern rankfun maxnum reverse
		      comparison-function lookupfun mapfun)
		     ;; not using up *
		     (trie--do-wildcard-search
		      node (trie--seq-append seq (trie--node-split node))
		      (cons ?* pattern) rankfun maxnum reverse
		      comparison-function lookupfun mapfun)))
		 (trie--node-subtree node)))

       ;; ? wildcard: map over all child nodes
       ((trie--wildcard-?-p el)
	(funcall mapfun
		 (lambda (node)
		   ;; skip data nodes (note: if we wanted to implement a "0 or
		   ;; 1" wildcard, would need to accumulate these instead)
		   (unless (trie--node-data-p node)
		     (trie--do-wildcard-search
		      node (trie--seq-append seq (trie--node-split node))
		      pattern rankfun maxnum reverse
		      comparison-function lookupfun mapfun)
		     ))
		 (trie--node-subtree node)
		 (if maxnum reverse (not reverse))))

       ;; character alternative: descend to corresponding nodes in turn
       ((trie--wildcard-char-alt-p el)
	(let (n)
	  (mapc
	   (lambda (c)
	     (when (setq n (funcall lookupfun (trie--node-subtree node)
				    (trie--node-create-dummy c)))
	       (trie--do-wildcard-search
		n (trie--seq-append seq c) pattern rankfun maxnum reverse
		comparison-function lookupfun mapfun)))
	   (if rankfun el
	     (sort el (if (or (and maxnum reverse)  ; no xnor in Elisp!
			      (and (not maxnum) (not reverse)))
			  (lambda (a b)
			    (not (funcall comparison-function a b)))
			comparison-function))))))

       ;; negated character alternative: map over all child nodes, skipping
       ;; excluded ones
       ((trie--wildcard-neg-char-alt-p el)
	(pop el)
	(funcall mapfun
		 (lambda (node)
		   ;; skip data nodes (note: if we wanted to implement a "0 or
		   ;; 1" wildcard, would need to accumulate these instead)
		   (unless (or (trie--node-data-p node)
			       (catch 'excluded
				 (dolist (c (butlast el))  ; drop final ^
				   (when (eq c (trie--node-split node))
				     (throw 'excluded t)))))
		     (trie--do-wildcard-search
		      node (trie--seq-append seq (trie--node-split node))
		      pattern rankfun maxnum reverse
		      comparison-function lookupfun mapfun)
		     ))
		 (trie--node-subtree node)
		 (if maxnum reverse (not reverse))))
       ))))



(defun trie-wildcard-stack (trie pattern &optional reverse)
  "blah"
  ;; convert trie from print-form if necessary
  (trie-transform-from-read-warn trie)
  ;; if stack functions aren't defined for trie type, throw error
  (if (not (functionp (trie--stack-createfun trie)))
      (error "Trie type does not support stack operations")
    ;; otherwise, create and initialise a stack
    (trie--wildcard-stack-create trie pattern reverse)))



(defun trie--wildcard-stack-construct-store
  (trie pattern &optional reverse)
  ;; Construct store for wildcard stack based on TRIE.
  (unless (or (atom pattern)
	      (and (listp pattern)
		   (not (sequencep (car pattern)))))
    (error "Multiple pattern searches are not currently supported by\
 trie-wildcard-stack's"))
  (let* ((comparison-function (trie--comparison-function trie))
	 (store
	 (list
	  (cons (cond ((stringp pattern) "") ((listp pattern) ()) (t []))
		;; convert pattern to list before parsing
		(cons
		 (trie--wildcard-parse-pattern
		  (append pattern nil)
		  (if reverse
		      `(lambda (a b) (,comparison-function b a))
		    comparison-function))
		 (trie--root trie))))))
    (trie--wildcard-stack-repopulate
     store reverse
     (trie--comparison-function trie)
     (trie--lookupfun trie)
     (trie--stack-createfun trie)
     (trie--stack-popfun trie)
     (trie--stack-emptyfun trie))))



(defun trie--wildcard-stack-repopulate
  (store reverse comparison-function lookupfun
	 stack-createfun stack-popfun stack-emptyfun)
  ;; Recursively push matching children of the node at the head of STORE onto
  ;; the front of STORE, until a data node is reached. Sort in (reverse)
  ;; lexical order if REVERSE is nil (non-nil). The remaining arguments should
  ;; be the corresponding trie functions (note that COMPARISON-FUNCTION should
  ;; be the trie--comparison-function, *not* the trie--cmpfun)
  (let (seq pattern node)
    (catch 'done
      (while t
	(catch 'cycle
	  ;; nothing to do if stack is empty
	  (unless store (throw 'done nil))


	  ;; if first stack element contains single node, and is not a character
	  ;; alternative, process it first
	  (setq seq (caar store)
		pattern (car (cdar store))
		node (cdr (cdar store)))
	  (when (trie--node-p node)
	    (setq store (cdr store))
	    ;; literal string: descend to corresponding node and continue
	    ;; processing (following element of pattern must be wildcard)
	    (when (trie--wildcard-literal-p (car pattern))
	      (setq node (trie--node-find node (car pattern) lookupfun))
	      ;; if we fail to find node corresponding to string, current
	      ;; branch of search has failed, so cycle and keep searching
	      (if (null node)
		  (throw 'cycle nil)
		;; if we found node corresponding to string, select that node
		(setq seq (trie--seq-concat seq (car pattern)))
		(setq pattern
		      (trie--wildcard-parse-pattern
		       (cdr pattern)
		       (if reverse
			   `(lambda (a b) (,comparison-function b a))
			 comparison-function)))))

	    (cond
	     ;; empty pattern: look for data node
	     ((null pattern)
	      (setq node (trie--find-data-node node lookupfun))
	      ;; if we fail to find one, current branch of search has failed,
	      ;; so cycle and keep searching
	      (if (null node)
		  (throw 'cycle nil)
		;; if we find one, push match onto stack and we're done
		(push (cons seq (trie--node-data node)) store)
		(throw 'done store)))

	     ;; character alternative: push node onto the stack
	     ((trie--wildcard-char-alt-p (car pattern))
	      (push (cons seq (cons pattern node)) store))

	     ;; any other wildcard: push a wildcard node stack onto the stack
	     (t (push (cons seq
			    (cons pattern
				  (funcall stack-createfun
					   (trie--node-subtree node) reverse)))
		      store))))


	  ;; first stack element is a wildcard pattern, so process it
	  (cond
	   ;; terminal *: standard repopulation using everything below node
	   ((and (null (cdr pattern)) (trie--wildcard-*-p (car pattern)))
	    ;; get first node from wildcard node stack
	    (setq node (funcall stack-popfun (cdr (cdar store))))
	    (when (funcall stack-emptyfun (cdr (cdar store)))
	      (setq store (cdr store)))
	    ;; recursively push node stacks for child nodes onto the stack until
	    ;; we find a data node
	    (while (not (trie--node-data-p node))
	      (push
	       (cons (trie--seq-append seq (trie--node-split node))
		     (cons pattern
			   (funcall stack-createfun
				    (trie--node-subtree node) reverse)))
	       store)
	      (setq node (funcall stack-popfun (cdr (cdar store)))
		    seq (caar store))
	      (when (funcall stack-emptyfun (cdr (cdar store)))
		(setq store (cdr store))))
	    (push (cons seq (trie--node-data node)) store)
	    (throw 'done store))

	   ;; non-terminal *: not currently supported
	   ((trie--wildcard-*-p (car pattern))
	    (error "Non-terminal * wildcards are not currently supported by\
 trie-wildcard-stack's"))

	   ;; ? wildcard: push wildcard node stack onto stack and repopulate
	   ;; again
	   ((trie--wildcard-?-p (car pattern))
	    ;; get first non-data node from wildcard node stack
	    (setq node (funcall stack-popfun (cdr (cdar store))))
	    (when (and node (trie--node-data-p node))
	      (setq node (funcall stack-popfun (cdr (cdar store)))))
	    (when (funcall stack-emptyfun (cdr (cdar store)))
	      (setq store (cdr store)))
	    (when node
	      (push
	       (cons (trie--seq-append seq (trie--node-split node))
		     (cons (trie--wildcard-parse-pattern
			    (cdr pattern)
			    (if reverse
				`(lambda (a b) (,comparison-function b a))
			      comparison-function))
			   node))
	       store)))

	   ;; character alternative: push next matching node onto stack and
	   ;; repopulate again
	   ((trie--wildcard-char-alt-p (car pattern))
	    (let ((c (pop (car pattern))))
	      (while (and c
			  (not (setq node
				     (funcall lookupfun
					      (trie--node-subtree node)
					      (trie--node-create-dummy c)))))
		(setq c (pop (car pattern))))
	      ;; if we've exhausted all characters in the alternative, remove it
	      ;; from the stack
	      (when (null (car pattern)) (setq store (cdr store)))
	      ;; if we found a match, push matching node onto stack and
	      ;; repopulate
	      (when node
		(push
		 (cons (trie--seq-append seq (trie--node-split node))
		       (cons (trie--wildcard-parse-pattern
			      (cdr pattern)
			      (if reverse
				  `(lambda (a b) (,comparison-function b a))
				comparison-function))
			     node))
		 store))))

	   ;; negated character alternative: push next non-excluded node onto
	   ;; stack and repopulate again
	   ((trie--wildcard-neg-char-alt-p (car pattern))
	    ;; pop nodes from wildcard node stack until we find one that isn't
	    ;; excluded
	    (setq node (funcall stack-popfun (cdr (cdar store))))
	    (while (and node
			(catch 'excluded
			  (dolist (c (butlast (car pattern))) ; drops final ^
			    (when (eq (trie--node-split node) c)
			      (throw 'excluded t)))))
	      (setq node (funcall stack-popfun (cdr (cdar store)))))
	    ;; remove wildcard node stack if empty
	    (when (funcall stack-emptyfun (cdr (cdar store)))
	      (setq store (cdr store)))
	    ;; if we found a match, push node onto stack; then repopulate again
	    (when node
	      (push
	       (cons (trie--seq-append seq (trie--node-split node))
		     (cons (trie--wildcard-parse-pattern
			    (cdr pattern)
			    (if reverse
				`(lambda (a b) (,comparison-function b a))
			      comparison-function))
			   node))
	       store)))
	   )
	  )))  ; end of infinite loop and catches
    )
  store)  ; return repopulated store



(defun trie--wildcard-parse-pattern (pattern &optional cmpfun)
  ;; Extract first pattern element from PATTERN (a list), and return it consed
  ;; with remainder of pattern. If CMPFUN is supplied, it is used to sort
  ;; character alternatives.
  (when pattern
    (let ((el (pop pattern)))
      (cond
       ;; *: drop any following *'s
       ((eq el ?*)
	(while (eq (car pattern) ?*) (pop pattern)))

       ;; [: gobble up to closing ]
       ((eq el ?\[)
	;; character alternatives are stored in lists
	(setq el ())
	;; gobble ] appearing straight after [
	(when (eq (car pattern) ?\]) (push (pop pattern) el))
	(while (not (eq (car pattern) ?\]))
	  (push (pop pattern) el)
	  (unless pattern
	    (error "Syntax error in trie wildcard pattern:\
 missing \"]\"")))
	(pop pattern)  ; dump closing ]
	;; if CMPFUN was supplied, sort characters in alternative
	(when cmpfun
	  ;; leave final ^ at end in negated character alternative
	  (if (eq (car (last el)) ?^)
	      (setq el (concat (sort (butlast el) cmpfun) ?^))
	    (setq el (sort el cmpfun)))))

       ;; ?: nothing to gobble
       ((eq el ??))

       ;; ]: syntax error (always gobbled when parsing [)
       ((eq el ?\])
	(error "Syntax error in trie wildcard pattern:\
 missing \"[\""))

       ;; anything else, gobble up to first special character
       (t
	(push el pattern)
	(setq el nil)
	(while (and pattern
		    (not (or (eq (car pattern) ?\[) (eq (car pattern) ?\])
			     (eq (car pattern) ?*) (eq (car pattern) ??))))
	  ;; \: dump \ and gobble next character
	  (when (eq (car pattern) ?\\)
	    (pop pattern)
	    (unless pattern
	      (error "Syntax error in trie wildcard pattern:\
 missing character after \"\\\"")))
	  (push (pop pattern) el))
	;; fixed strings are stored in vectors
	(setq el (vconcat (nreverse el)))))

      ;; return cons containing first element and remaining pattern
      (cons el pattern))))



(provide 'trie)

;;; trie.el ends here
