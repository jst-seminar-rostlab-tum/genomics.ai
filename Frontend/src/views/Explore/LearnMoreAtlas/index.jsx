import React, { useEffect, useState } from 'react';
import NavBar from 'components/NavBar';
import { Box, Button, Typography } from '@mui/material';
import Breadcrumb from 'components/Breadcrumb';
import Chip from '@mui/material/Chip';
import CustomButton from 'components/CustomButton';
import Search from 'components/Search';
import { Filter } from 'components/Filter/Filter';
import AtlasService from 'shared/services/Atlas.service';
import { useLocation } from 'react-router-dom';

const mockData = {
  name: 'Human - PBMC', previewPictureURL: 'url', modalities: ['RNA', 'ADT'], numberOfCells: 161.764, species: ['Human'], compatibleModels: ['ObjectId'],
};
// FIXME: Those are temproary, Tag and DataSet components should be fix then
// replace them with temproaryDataSetCards below.
const TemproaryDataSetCard = () => (
  <div style={{
    height: '41px',
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    background: 'white',
    borderRadius: '10px',
    fontSize: '13px',
    width: '300px',
    heigth: '41px',
    border: '1px solid white',
    filter: 'drop-shadow(0px 0px 7px rgba(0, 0, 0, 0.25))',
  }}
  >
    <div style={{
      display: 'flex', alignItems: 'center', columnGap: '8px', height: '21px',
    }}
    >
      <p>Hao and Hao et al, bioRvix 2020</p>
      <span style={{
        display: 'flex',
        borderRadius: '10px',
        fontSize: '13px',
        height: '21px',
        color: 'white',
        background: '#01579B',
        alignItems: 'center',
        padding: '4px',
      }}
      >
        category
      </span>

    </div>
  </div>
);

const TemproaryDataSetCard2 = () => (
  <div style={{
    height: '41px',
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    background: 'white',
    borderRadius: '10px',
    fontSize: '13px',
    width: '300px',
    heigth: '41px',
    border: '1px solid white',
    filter: 'drop-shadow(0px 0px 7px rgba(0, 0, 0, 0.25))',
  }}
  >
    <div style={{
      display: 'flex', alignItems: 'center', columnGap: '8px', height: '21px',
    }}
    >
      <p>11,769 PBMCs from 10x Genomics</p>
      <span style={{
        display: 'flex',
        borderRadius: '10px',
        fontSize: '13px',
        height: '21px',
        color: 'white',
        background: '#01579B',
        alignItems: 'center',
        padding: '4px',
      }}
      >
        category
      </span>

    </div>
  </div>
);

export default function LearnMore() {
  const [value, setValue] = useState(0);
  const path = useLocation().pathname;
  const id = path.match("/explore/atlases/(.*)$")[1]
  const [atlas, setAtlas] = useState(null);

  /* const [selectedAtlas, setSelectedAtlas] = useState(null);
  const [selectedModel, setSelectedModel] = useState(null); */
  const [mapperVisible, setMapperVisible] = useState(false);

  useEffect(() => {
    if(id) {
      AtlasService.getAtlasById(id)
        .then((data) => setAtlas(data))
        .catch((err) => console.log(err));
    }
  }, [id]);


  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>

      <Box sx={{
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'flex-start',
        marginTop: '12px',
        width: '80%',
        minWidth: '1200px',
        justifyContent: 'space-between',
      }}
      >
        <Box sx={{
          display: 'flex', flexDirection: 'row', width: '100%', justifyContent: 'space-between',
        }}
        >
          <Typography sx={{ fontSize: '36px', fontWeigth: 700 }}>{atlas?.name}</Typography>
        </Box>
        <Box>
          <Typography sx={{ fontSize: '20px', fontWeight: 600, borderBottom: '1px solid black' }}>Overview</Typography>
        </Box>
        <Box sx={{ display: 'flex', flexDirection: 'row', paddingTop: '16px' }}>
          <Typography sx={{ fontSize: '16px', fontWeight: 500 }}>
            Modalities:
            &nbsp;
          </Typography>
          <Typography sx={{ fontSize: '16px', fontWeight: 300 }}>
            {atlas?.modalities}
          </Typography>
        </Box>
        <Box sx={{ display: 'flex', flexDirection: 'row' }}>
          <Typography sx={{ fontSize: '16px', fontWeight: 500 }}>
            Cells in Reference:
            &nbsp;
          </Typography>
          <Typography sx={{ fontSize: '16px', fontWeight: 300 }}>
            {atlas?.numberOfCells}
          </Typography>
        </Box>
        <Box sx={{ display: 'flex', flexDirection: 'row' }}>
          <Typography sx={{ fontSize: '16px', fontWeight: 500 }}>
            Species:
            &nbsp;
          </Typography>
          <Typography sx={{ fontSize: '16px', fontWeight: 300 }}>
            {atlas?.species}
          </Typography>
        </Box>
        
        {/* <Box sx={{ display: 'flex', columnGap: '18px', paddingTop: '18px' }}>
          <Chip label="App" sx={{ background: '#01579B', color: 'white', height: '26px' }} />
          <Chip label="Reference" sx={{ background: '#01579B', color: 'white', height: '26px' }} />
          <Chip label="Zenodo" sx={{ background: '#01579B', color: 'white', height: '26px' }} />
          <Chip label="Snakemake" sx={{ background: '#01579B', color: 'white', height: '26px' }} />
        </Box> */}
        {/* <Box sx={{
          display: 'flex', flexDirection: 'column', width: '300px', rowGap: '18px', paddingTop: '18px',
        }}
        >
          <Typography sx={{
            fontSize: '20px', textDecoration: 'underline solid black 1px', textUnderlineOffset: 10, fontWeight: 600,
          }}
          >
            Reference Dataset(s)
          </Typography>
          <TemproaryDataSetCard />
        </Box>
        <Box sx={{
          display: 'flex', flexDirection: 'column', width: '300px', rowGap: '18px', paddingTop: '18px',
        }}
        >
          <Typography sx={{
            fontSize: '20px', textDecoration: 'underline solid black 1px', textUnderlineOffset: 10, fontWeight: 600,
          }}
          >
            Demo Dataset(s)
          </Typography>
          <TemproaryDataSetCard2 />
        </Box>
        <Box sx={{
          display: 'flex', flexDirection: 'column', width: '300px', rowGap: '18px', paddingTop: '18px',
        }}
        >
          <Typography sx={{
            fontSize: '20px', fontWeight: 600, textDecoration: 'underline solid black 1px', textUnderlineOffset: 10,
          }}
          >
            Annotation Details
          </Typography>
        </Box> */}
        <CustomButton sx={{ marginTop: "1em", padding: "0.5em 2em 0.5em 2em" }} type="primary">Map</CustomButton>
      </Box>

    </Box>
  );
}
