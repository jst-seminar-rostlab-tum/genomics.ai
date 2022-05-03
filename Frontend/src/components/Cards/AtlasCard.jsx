/* eslint-disable react/react-in-jsx-scope */
import { useState, useRef, useEffect } from 'react';

import { Box, Typography } from '@mui/material';

import { OutlinedButton } from './ModelCard';
import { useHistory } from 'react-router-dom/cjs/react-router-dom.min';

/**
 * Atlas Card
 * @param width default value is 100% of parent
 * @param height default value is 100% of parent
 * @param title title of AtlasCard
 * @param imgLink thumbnail photo url
 * @param modalities
 * @param cellsInReference
 * @param species
 * @param mapLink onHover button Map url
 * @param learnMoreLink onHover button Learn More url
 */
export default function AtlasCard({
  width = "100%", height = "100%", title, atlasId, imgLink, modalities,
  cellsInReference, species, learnMoreLink, onClick, atlas
}) {
  // check if the mouse is hovering above the card
  const [isHover, setHover] = useState(false);

  // check if the card is flat(width > height)
  const [isFlat, setFlat] = useState(false);

  // ref to get the out most Box
  const boxRef = useRef();

  const history = useHistory();
  const path = history.location.pathname;

  useEffect(() => {
    //each time the card is rerendered, check if the card is flat or not
    if (boxRef.current.clientWidth > boxRef.current.clientHeight) setFlat(true)

  })

  const showModalities = () => {
    //TODO fix this
    if (modalities[0].length < 10) return modalities
    return `${modalities[0].split(",")[0]}`
  }

  const AtlasInfo = (title, data) => (
    <Box
      sx={{
        display: "flex",
        flexDirection: data.length > 10 ? "column" : 'row',
        gap: "5%"
      }}
    >
      <Typography
        sx={{
          fontSize: "1rem",
          fontWeight: "bold"
        }}
      >
        {title}
      </Typography>
      <Typography>{data}</Typography>
    </Box>
  )

  return (
    <Box
      sx={{
        width, height,
      }}
    >
      <Box
        ref={boxRef}
        onMouseEnter={() => setHover(true)}
        onMouseLeave={() => setHover(false)}
        sx={{
          width: '100%',
          height: '100%',
          position: 'relative',
        }}
      >
        {
          isHover
          && (
            <Box
              style={{
                background: 'linear-gradient(#4F83CC, #01579B)',
              }}
              sx={{
                position: "absolute",
                width: "100%",
                height: "100%",
                borderRadius: "1.2rem",
                display: "flex",
                flexDirection: "column",
                justifyContent: "center",
                opacity: 0.95,
                boxShadow: "0px 4px 6px 0px rgba(1, 87, 155, .20), 0px 0px 1px 0px rgba(1, 87, 155, .32)"
              }}
            >
              <Box
                sx={{
                  margin: 'auto',
                  width: isFlat ? '70%' : '70%',
                  height: isFlat ? 'auto' : '40%',
                  display: 'flex',
                  flexDirection: isFlat ? 'row' : 'column',
                  justifyContent: 'space-evenly',
                }}
              >
                <OutlinedButton content="Map" onClick={onClick} />
                <OutlinedButton content="Learn More" link={learnMoreLink} onClick={() => localStorage.setItem("atlasId", atlasId)} />
              </Box>
            </Box>
          )
        }

        <Box
          sx={{
            width: "100%",
            height: "100%",
            padding: "1rem",
            display: "flex",
            flexDirection: "column",
            boxShadow: isHover ? 'none' : "0px 4px 6px 0px rgba(33, 37, 41, .2), 0px 0px 1px 0px rgba(33, 37, 41, .32)",
            borderRadius: "1.2rem"
          }}
        >
          <Typography
            sx={{
              fontSize: '1.4rem',
              fontWeight: 'bold',
            }}
          >
            {title}
          </Typography>

          <Box
            component="img"
            src={imgLink}
            alt="Atlas preview img"
            sx={{
              width: "90%",
              objectFit: "cover", 
              alignSelf: 'center'
            }}
          />
          <Box sx={{ display: 'flex', flexDirection: 'column', m: "5px ", justifyContent:'space-evenly', height:'100%'}}>
            {AtlasInfo("Modalities", showModalities())}
            {AtlasInfo("Species", species)}
          </Box>

        </Box>
      </Box>
    </Box>
  );
}
