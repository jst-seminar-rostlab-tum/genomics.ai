import { useState, useRef, useEffect } from "react"

import { Box, Container, Typography } from '@mui/material'

import { OutlinedButtonSelect } from './ModelCardSelect'

import { borders } from "@mui/system"
import { Modal } from "components/Modal"
import { LearnMoreAtlasComponent } from "views/Explore/LearnMoreAtlas"
import { useHistory } from "react-router-dom/cjs/react-router-dom.min"
import AtlasInfo from "components/GeneMapper/AtlasInfo"

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
export default function AtlasCardSelect({
  width = "80", height = "80%", title, imgLink, modalities,
  cellsInReference, species, mapLink, learnMoreLink, selected=false, 
  onSelect, atlasObject={}
}) {

  //check if the mouse is hovering above the card
  const [isHover, setHover] = useState(false)

  //check if the card is flat(width > height)
  const [isFlat, setFlat] = useState(false)


  const [atlasInfoOpen, setAtlasInfoOpen] = useState(false);

  // ref to get the out most Box
  const boxRef = useRef()
  const history = useHistory();

  useEffect(() => {
    // each time the card is rerendered, check if the card is flat or not
    if (boxRef.current.clientWidth > boxRef.current.clientHeight) setFlat(true)
  }, []);

  // capitalize every word in a string
  const capitalize = (str) => {
    const lower = str.toLowerCase();
    const words = lower.split(' ');
    // return the original string if it was max one word since it is already capitalized and correct
    console.log(words);
    if (words.length === 1) return str;
    let cap = '';
    // eslint-disable-next-line no-plusplus
    for (let i = 0; i < words.length; i++) {
      const word = words[i].toLowerCase();
      // split the word into the first letter and the rest of the string
      let firstLetter = word.charAt(0);
      firstLetter = firstLetter.toUpperCase();
      const restOfWord = word.slice(1);
      cap += firstLetter + restOfWord;
      // add a space after every word except the last one
      if (i < words.length - 1) cap += ' ';
    }
    return cap;
  };

  return (
    <Box
      sx={{
        width, height
      }}
    >
      <Box
        ref={boxRef}
        onMouseEnter={() => setHover(true)}
        onMouseLeave={() => setHover(false)}
        sx={{
          width: "80%",
          height: "80%",
          position: "relative",
        }}
      >
        {/*Hover effect over the card*/}
        {
          isHover &&
          <Box
            style={{
              background: "linear-gradient(#4F83CC, #01579B)"
            }}
            sx={{
              position: "absolute",
              width: "100%",
              height: "100%",
              borderRadius: "1.2rem",
              display: "flex",
              flexDirection: "column",
              justifyContent: "center",
              opacity: 0.95
            }}
          >
            <Box
              sx={{
                margin: "auto",
                width: isFlat ? "70%" : "60%",
                height: isFlat ? "auto" : "40%",
                display: "flex",
                flexDirection: isFlat ? "row" : "column",
                justifyContent: "space-evenly",
              }}
            >
              <OutlinedButtonSelect content="Select" onSelect={() => onSelect(atlasObject)} />
              <OutlinedButtonSelect content="Learn More" onSelect={() => setAtlasInfoOpen(true)} />
            </Box>
          </Box>
        }


        <Box
          sx={{
            width: '100%',
            height: '100%',
            padding: '1rem',
            display: 'flex',
            flexDirection: 'column',
            boxShadow: isHover ? 'none' : '0px 4px 6px 0px rgba(33, 37, 41, .2), 0px 0px 1px 0px rgba(33, 37, 41, .32)',
            borderRadius: '1.2rem',
            justifyContent: 'center',
            borderStyle: 'solid',
            borderColor: selected ? '#008BF5' : 'transparent',
            borderWidth: '4px',
          }}
        >
          <Typography
            sx={{
              fontSize: "1.2rem",
              fontWeight: "bold"
            }}
          >
            {capitalize(title)}
          </Typography>

          <Box component="img" src={imgLink} alt="Atlas preview img"
            sx={{
              width: "90%",
              height: "50%",
              margin: "auto"
            }}
          />

          <Box
            sx={{
              display: "flex",
              flexDirection: "row",
            }}
          >
            <Typography
              sx={{
                fontSize: "1rem",
                fontWeight: "bold"
              }}
            >
              Modalities:
            </Typography>
            &nbsp;
            <Typography noWrap>{modalities}</Typography>
          </Box>

          <Box
            sx={{
              display: "flex",
              flexDirection: "row",
            }}
          >
            <Typography
              sx={{
                fontSize: "1rem",
                fontWeight: "bold"
              }}
            >
              Cells in Reference:
            </Typography>
            &nbsp;
            <Typography>{cellsInReference}</Typography>
          </Box>

          <Box
            sx={{
              display: "flex",
              flexDirection: "row",
            }}
          >
            <Typography
              sx={{
                fontSize: "1rem",
                fontWeight: "bold"
              }}
            >
              Species:
            </Typography>
            &nbsp;
            <Typography>{species}</Typography>
          </Box>
        </Box>
      </Box>
      <AtlasInfo id={atlasObject._id} open={atlasInfoOpen} setOpen={setAtlasInfoOpen} />
    </Box>

  )
}