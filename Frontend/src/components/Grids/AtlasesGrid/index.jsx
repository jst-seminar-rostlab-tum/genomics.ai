import React from 'react';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import AtlasCard from 'components/Cards/AtlasCard';
import styles from './atlasesGrid.module.css';

const AtlasesGrid = ({
  atlases, path, selectedAtlas, setSelectedAtlas,
}) => (
  <Box className={styles.atlasContainer} maxHeight="50vh" mb="2em">
    <Grid container spacing={3}>
      {atlases && atlases.map((atlas) => (
        <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
          <AtlasCard
            atlasId={atlas._id}
            imgLink={atlas.previewPictureURL}
            species={atlas.species}
            modalities={atlas.modalities}
            title={atlas.name}
            learnMoreLink={`${path}/atlases/${atlas._id}`}
            onSelect={() => setSelectedAtlas(atlas)}
            selected={selectedAtlas && selectedAtlas.name === atlas.name}
          />
        </Grid>
      ))}
    </Grid>
  </Box>
);

export default AtlasesGrid;
