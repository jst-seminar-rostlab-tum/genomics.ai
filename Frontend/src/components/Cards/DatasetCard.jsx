/* eslint-disable */

import { useEffect, useRef, useState } from 'react'

import { Box, Typography } from '@mui/material'

import Tag from 'components/Tag'

/**
 * Dataset Card
 * @param width default value 100%
 * @param height default value 100%
 * @param title
 * @param category
 */
export default function DatasetCard({ title, category, width = "100%", height = "100%" }) {

  //store the size of the card in order to set the size of tag properly
  const [cardSize, setCardSize] = useState({ width: 0, height: 0 })

  //store the ref of the card
  const cardRef = useRef()

  useEffect(() => {
    //before rendered, store the size of the card
    setCardSize({ width: cardRef.current.clientWidth, height: cardRef.current.clientHeight })
  }, [])

  return (
    <Box
      sx={{
        width: width,
        height: height
      }}
    >
      <Box ref={cardRef}
        sx={{
          width: "100%",
          height: "100%",
          display: "flex",
          flexDirection: "row",
          justifyContent: "space-evenly",
          padding: "5px",
          boxShadow: "0px 4px 6px 0px rgba(33, 37, 41, .2), 0px 0px 1px 0px rgba(33, 37, 41, .32)",
          borderRadius: "10px"
        }}
      >
        <Box
          sx={{
            display: "flex",
            flexDirection: "column",
            justifyContent: "center"
          }}
        >
          <Typography fontSize=".9rem" fontWeight="bold">
            {title}
          </Typography>
        </Box>

        <Box
          sx={{
            display: "flex",
            flexDirection: "column",
            justifyContent: "center"
          }}
        >
          <Box
            sx={{
              width: `${0.25 * cardSize.width}px`,
              height: `${0.45 * cardSize.height}px`
            }}
          >
            <Tag content={category} fontSize={13} variant="secondary-default" isDeletable={false} />
          </Box>
        </Box>
      </Box>
    </Box>
  )
}
