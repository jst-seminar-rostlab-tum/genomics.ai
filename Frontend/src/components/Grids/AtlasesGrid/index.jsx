import React from 'react';
import { useHistory, useLocation } from 'react-router-dom';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import AtlasCard from 'components/Cards/AtlasCard';
import styles from './atlasesGrid.module.css';

function applyAtlasFilters(atlases, searchedKeyword, searchParams) {
  const searchedAtlases = atlases.filter(
    (item) => item.name.toLowerCase().includes(searchedKeyword.toLowerCase()),
  );
  if (searchParams.get('sortBy') === 'name' || searchParams.get('sortBy') === null) {
    searchedAtlases.sort((a, b) => {
      const nameA = a.name.toUpperCase();
      const nameB = b.name.toUpperCase();
      if (nameA < nameB) {
        return -1;
      }
      if (nameA > nameB) {
        return 1;
      }
      return 0;
    });
  } else if (searchParams.get('sortBy') === 'numberOfCells') {
    searchedAtlases.sort((a, b) => a.numberOfCells - b.numberOfCells);
  }
  return searchedAtlases;
}

const AtlasesGrid = ({ atlases, searchedKeyword, path }) => {
  const history = useHistory();
  const { search } = useLocation();
  const searchParams = new URLSearchParams(search);
  const atlasesFiltered = applyAtlasFilters(atlases, searchedKeyword, searchParams);
  return (
    <Box className={styles.atlasContainer} sx={{ height: '70vh' }}>
      <Grid container spacing={3}>
        {atlasesFiltered && atlasesFiltered.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => history.push(`${path}/atlases/${atlas._id}/visualization`)}
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/atlases/${atlas._id}`}
            />
          </Grid>
        ))}
        {atlasesFiltered && atlasesFiltered.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => history.push(`${path}/atlases/${atlas._id}/visualization`)} // removed setSelectedAtlas=> feels outdated
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/atlases/${atlas._id}`}
            />
          </Grid>
        ))}
        {atlasesFiltered && atlasesFiltered.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => history.push(`${path}/atlases/${atlas._id}/visualization`)} // removed setSelectedAtlas=> feels outdated
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/atlases/${atlas._id}`}
            />
          </Grid>
        ))}
      </Grid>
    </Box>
  );
};

export default AtlasesGrid;
