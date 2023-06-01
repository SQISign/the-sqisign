;; Emacs style file
;;
;; Sets spaces-only indentation, 4-spaces tab stops, linux kernel
;; coding style
(
 (nil . ((indent-tabs-mode . nil)
         (tab-width . 4)
         )
      )
 (c-default-style . ((c-mode . "linux")
                     ))
 (c-mode . ((c-file-style . "linux")
            (c-basic-offset . 4)
            )
         )
 )
