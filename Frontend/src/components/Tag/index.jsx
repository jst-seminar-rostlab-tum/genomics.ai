/* eslint-disable */

import { useEffect, useRef, useState } from 'react'

import { Box, Typography } from '@mui/material'
import { Cancel } from '@mui/icons-material'

import { colors } from 'shared/theme/colors'

/**
 * predefined Tags for variant
 * 
 * if you want to customize a tag, you can use a object as variant and must contains the attribute:
 * - backgroundColor
 * - borderColor
 * - borderWidth
 * - fontColor
 */
const predefinedTags = [
  {
    variant: "primary-default",
    backgroundColor: "primary.main",
    borderColor: "primary.main",
    borderWidth: 0,
    fontColor: colors.neutral[100]
  },
  {
    variant: "secondary-default",
    backgroundColor: "secondary.main",
    borderColor: "secondary.main",
    borderWidth: 0,
    fontColor: colors.neutral[100]
  },
  {
    variant: "neutral-outlined",
    backgroundColor: "white",
    borderColor: colors.neutral[700],
    borderWidth: 2,
    fontColor: colors.neutral[700]
  }
]

/**
 * Tag
 * needs content, variant(see above), whether it's deletable, meaning whether the delete icon is to shown and a function to handle the deletion
 */
export default function Tag({ content, variant, fontSize=16, isDeletable, handleDeletion = () => { } }) {

  //store the offset of the delete icon to make its center right on the corner of Tag
  const [deleteIconOffset, setDeleteIconOffset] = useState({ x: 0, y: 0 })

  //store the ref to the delete icon
  const deleteIconRef = useRef()

  //before rendered, get the size of the delete button and set the offset correspondingly
  useEffect(() => {
    if (isDeletable) setDeleteIconOffset({ x: deleteIconRef.current.clientWidth / 3, y: -deleteIconRef.current.clientHeight / 3 })
  }, [])

  /**
   * the 4 get functions
   * 
   * if variant is a string, then it represents a predefined tag, thus search in the array
   * if variant is a object, get the corresponding attribute from it
   * otherwise do not do anything
   */
  function getBackgroundColor() {
    if (typeof variant === "string") {
      for (let tag of predefinedTags) {
        if (tag.variant === variant) return tag.backgroundColor
      }
      return predefinedTags[0].backgroundColor
    } else if (typeof variant === "object")
      return variant.backgroundColor
  }

  function getBorderColor() {
    if (typeof variant === "string") {
      for (let tag of predefinedTags) {
        if (tag.variant === variant) return tag.borderColor
      }
      return predefinedTags[0].borderColor
    } else if (typeof variant === "object")
      return variant.borderColor
  }

  function getBorderWidth() {
    if (typeof variant === "string") {
      for (let tag of predefinedTags) {
        if (tag.variant === variant) return tag.borderWidth
      }
      return predefinedTags[0].borderWidth
    } else if (typeof variant === "object")
      return variant.borderWidth
  }

  function getFontColor() {
    if (typeof variant === "string") {
      for (let tag of predefinedTags) {
        if (tag.variant === variant) return tag.fontColor
      }
      return predefinedTags[0].fontColor
    } else if (typeof variant === "object")
      return variant.fontColor
  }

  return (
    <Box
      sx={{
        position: "relative",
        width: "100%",
        height: "100%",
        border: getBorderWidth(),
        borderColor: getBorderColor(),
        backgroundColor: getBackgroundColor(),
        borderRadius: "20px",
      }}
    >
      {
        isDeletable &&
        <Box onClick={handleDeletion}
          sx={{
            position: "absolute",
            display: "flex",
            flexDirection: "row-reverse",
            width: "100%",
            height: "100%",
            left: deleteIconOffset.x,
            top: deleteIconOffset.y,
          }}
        >
          <Cancel ref={deleteIconRef}
            sx={{
              color: colors.error.main,
              bgcolor: "white",
              borderRadius: "100%"
            }}
          >
            {
              isDeletable &&
              <Box onClick={handleDeletion}
                sx={{
                  position: "absolute",
                  display: "flex",
                  flexDirection: "row-reverse",
                  width: "100%",
                  height: "100%",
                  left: deleteIconOffset.x,
                  top: deleteIconOffset.y,
                }}
              >
                <Cancel ref={deleteIconRef}
                  sx={{
                    color: colors.error.main,
                    bgcolor: "white",
                    borderRadius: "100%"
                  }}
                />
              </Box>
            }

            <Box
              sx={{
                display: "flex",
                flexDirection: "column",
                justifyContent: "center",
                width: "100%",
                height: "100%"
              }}
            >
              <Typography fontSize="70%" fontWeight="bold"
                sx={{
                  color: getFontColor(),
                  margin: "auto"
                }}
              >
                {content}
              </Typography>
            </Box>
          </Cancel>
        </Box>
      }

      <Box
        sx={{
          display: "flex",
          flexDirection: "column",
          justifyContent: "center",
          width: "100%",
          height: "100%"
        }}
      >
        <Typography fontSize="100%" fontWeight="bold"
          sx={{
            color: getFontColor(),
            margin: "auto",
            fontSize:{fontSize}
          }}
        >
          {content}
        </Typography>
      </Box>
    </Box>
  )
}
