import React from 'react';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import AtlasCard from 'components/Cards/AtlasCard';
import styles from './atlasesGrid.module.css';

const AtlasesGrid = ({
  atlases, path, selectedAtlas = null, handleAtlasSelection = null, selectedModel = null, isSearchPage = false
}) => (
  <Box className={styles.atlasContainer} mb="2em" pb="1.5em"
    sx={{
      
    }}
  >
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
            onSelect={() => { if (handleAtlasSelection) handleAtlasSelection(atlas); }}
            selected={selectedAtlas && selectedAtlas.name === atlas.name}
            disabled={selectedModel && !atlas.compatibleModels.some(
              (element) => element.toLowerCase() === selectedModel.name.toLowerCase(),
            )}
            isSearchPage={isSearchPage}
          />
        </Grid>
      ))}
    </Grid>
  </Box>
);

export default AtlasesGrid;
