/* eslint-disable react/react-in-jsx-scope */
import { useState, useRef, useEffect } from 'react';
import { Box, Typography } from '@mui/material';
import { OutlinedButton } from './ModelCard';
import { colors } from 'shared/theme/colors';
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
  width = '100%', height = '100%', title, atlasId, imgLink, modalities,
  cellsInReference, species, learnMoreLink, onSelect, selected = false, disabled = false,
  isSearchPage = false
}) {
  // check if the mouse is hovering above the card
  const [isHover, setHover] = useState(false);

  // check if the card is flat(width > height)
  {/*const [isFlat, setFlat] = useState(false);*/}

  // ref to get the out most Box
  const boxRef = useRef();
  const history = useHistory();
  const path = history.location.pathname;
  {/*useEffect(() => {
    // each time the card is rerendered, check if the card is flat or not
    if (boxRef.current.clientWidth > boxRef.current.clientHeight) setFlat(true);
  });*/}

  const showModalities = () => {
    // TODO fix this
    if (modalities[0].length < 10) return modalities;
    return `${modalities[0].split(',')[0]}`;
  };

  const AtlasInfo = (title, data) => (
    <Box
      sx={{
        display: 'flex',
        flexDirection: data.length > 10 ? 'column' : 'row',
        gap: '5%',
      }}
    >
      <Typography
        sx={{
          fontSize: '1rem',
          fontWeight: 'bold',
        }}
      >
        {title}
      </Typography>
      <Typography>{data}</Typography>
    </Box>
  );

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
          cursor: disabled ? 'default' : 'pointer',
        }}
      >
        {
          !disabled && isHover
          && (
            <Box
              style={{
                background: 'linear-gradient(#4F83CC, #01579B)',
              }}
              sx={{
                position: 'absolute',
                width: '100%',
                height: '100%',
                borderRadius: '1.2rem',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'center',
                cursor: 'pointer',
                opacity: 0.95,
                boxShadow: '0px 4px 6px 0px rgba(1, 87, 155, .20), 0px 0px 1px 0px rgba(1, 87, 155, .32)',
              }}
            >
              <Box
                sx={{
                  margin: 'auto',
                  width: '70%',
                  height: 'auto',
                  display: 'flex',
                  flexDirection: 'column',
                  justifyContent: 'space-evenly',
                  gap: '5px',
                }}
              >
                {
                  !isSearchPage &&
                  <OutlinedButton
                    content="Select"
                    onClick={onSelect}
                  />
                }
                <OutlinedButton
                  content="Learn More"
                  link={learnMoreLink}
                  onClick={() => {
                    localStorage.setItem('atlasId', atlasId);
                  }}
                />
                {
                disabled
                && (
                <Typography sx={{
                  color: colors.primary[900], fontSize: '12px', textDecoration: 'underline', textAlign: 'center',
                }}
                >
                  Atlas unavailable for mapping with selected model
                </Typography>
                )
              }
              </Box>
            </Box>
          )
        }

        {/* DISABLED OVERLAY BOX */}
        {
          disabled &&
          <Box
            style={{ background: "linear-gradient(#e7e7e7, #d0d0d0)" }}
            sx={{
              position: "absolute",
              width: "100%",
              height: "100%",
              display: 'flex',
              flexDirection: 'column',
              justifyContent: 'center',
              borderRadius: "1.2rem",
              cursor: "pointer",
              opacity: 0.95,
            }}>
            { 
              !isHover && 
              <Typography sx={{ position: "absolute", fontSize: "12px", fontWeight: "bold", color: colors.neutral[900], textAlign: "center", left: "28%" }}>
                Not Compatible
              </Typography>
            }
            <Box
              sx={{
                margin: 'auto',
                width: '70%',
                height: 'auto',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-evenly',
                gap: '5px',
              }}
            >

              { isHover && <OutlinedButton content="Learn More" link={learnMoreLink} bg={colors.neutral[800]} color={colors.neutral[800]} onClick={(e) => e.stopPropagation()} bgHover={colors.neutral[100]} colorHover={colors.neutral[600]}/> }
              {/* <CustomButton type="outlined" href={learnMoreLink ? `#${learnMoreLink}` : null} onClick={(e) => e.stopPropagation()}>Learn More</CustomButton> */}
            </Box>
          </Box>
        }
        { !disabled
          && (
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
                width: '90%',
                objectFit: 'cover',
                alignSelf: 'center',
              }}
            />
            <Box sx={{
              display: 'flex', flexDirection: 'column', m: '5px ', justifyContent: 'space-evenly', height: '100%',
            }}
            >
              {AtlasInfo('Modalities', showModalities())}
              {AtlasInfo('Species', species)}
            </Box>
          </Box>
          )}
        {
          disabled
          && (
          <Box sx={{
            width: '100%',
            height: '100%',
            display: 'flex',
            flexDirection: 'column',
            p: '1.2rem',
            borderRadius: '1.2rem',
          }}
          >
            <Typography sx={{ fontSize: '1.4rem', fontWeight: 'bold' }}>{title}</Typography>
            <Box
              component="img"
              src={imgLink}
              alt="Atlas preview img"
              sx={{
                width: '90%',
                objectFit: 'cover',
                alignSelf: 'center',
              }}
            />
            <Box sx={{
              display: 'flex', flexDirection: 'column', m: '5px ', justifyContent: 'space-evenly', height: '100%',
            }}
            >
              {AtlasInfo('Modalities', showModalities())}
              {AtlasInfo('Species', species)}
            </Box>
          </Box>
          )
        }
      </Box>
    </Box>
  );
}
