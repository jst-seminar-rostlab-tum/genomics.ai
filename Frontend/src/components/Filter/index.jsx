import React, { useState } from 'react';
import FormGroup from '@mui/material/FormGroup';
import Box from '@mui/material/Box';
import MenuItem from '@mui/material/MenuItem';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import FormControlLabel from '@mui/material/FormControlLabel';
import Checkbox from '@mui/material/Checkbox';
import styles from './filter.module.css';
import { Divider, Stack, TextField } from '@mui/material';
import Tag from '../Tag';

export default function Filter() {
  const [reference, setReference] = useState('');
  const [category, setCategory] = useState('');
  const [selectedIndexes, setSelectedIndexes] = useState([]);
  const references = ['Human - PBMC', 'Human - PBMC2', 'Human - PBMC3', 'Human - PBMC4', 'Human - PBMC5'];
  const categories = ['lung', 'lung2', 'lung3', 'lung4', 'lung5'];

  const handleReferenceChange = (event) => {
    setReference(event.target.value);
  };

  const handleCategoryChange = (event) => {
    setCategory(event.target.value);
  };

  const handleIndexesAppend = (newValue) => {
    setSelectedIndexes((array) => [...array, newValue]);
  };

  const addCategory = () => {
    if (!selectedIndexes.includes(category) && category !== '') {
      handleIndexesAppend(category);
    }
  };

  const removeTag = (selectedIndex) => {
    setSelectedIndexes(selectedIndexes.filter((item) => item !== selectedIndex));
  };

  return (
    <Box className={styles.container}>
      <FormGroup>
        <Stack
          direction="row"
          className={styles.stack}
        >
          <FormControlLabel
            control={<Checkbox />}
            label={(
              <Typography
                className={styles.checkboxLabel}
              >
                Referenced By Atlas
              </Typography>
)}
          />
          <TextField
            autoWidth
            value={reference}
            onChange={handleReferenceChange}
            className={styles.select}
            select
            variant="standard"
          >
            {references.map((referenceItem, index) => (
              <MenuItem value={index}>{referenceItem}</MenuItem>
            ))}
          </TextField>
        </Stack>
        <Divider className={styles.divider} />
        <Stack
          direction="row"
          className={styles.stack}
        >
          <FormControlLabel
            control={<Checkbox />}
            label={(
              <Typography
                className={styles.checkboxLabel}
              >
                Category
              </Typography>
)}
          />
          <TextField
            autoWidth
            value={category}
            onChange={handleCategoryChange}
            className={styles.select}
            select
            variant="standard"
          >
            {categories.map((categoryItem, index) => (
              <MenuItem value={index}>{categoryItem}</MenuItem>
            ))}
          </TextField>
          <Button variant="contained" className={styles.addButton} onClick={addCategory}>Add</Button>
        </Stack>
        <Stack
          direction="row"
          className={styles.tagStack}
        >
          {selectedIndexes.map((selectedIndex) => (
            <Tag text={categories[selectedIndex]} handleClick={() => removeTag(selectedIndex)} />
          ))}
        </Stack>
      </FormGroup>

    </Box>
  );
}
