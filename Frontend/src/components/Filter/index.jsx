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

export default function Filter(props) {
  const [reference, setReference] = useState('');
  const [category, setCategory] = useState('');
  const [selectedIndexes, setSelectedIndexes] = useState([]);
  const { references, categories } = props;

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
    <Box className={styles.container}
      sx={{
        backgroundColor: "white",
        borderRadius: "20px",
        boxShadow: "0px 4px 6px rgba(0, 0, 0, 0.15), 0px 0px 1px rgba(0, 0, 0, 0.4)",
        zIndex: "100"
      }}
    >
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
            value={reference}
            onChange={handleReferenceChange}
            className={styles.select}
            select
            variant="standard"
          >
            {references.map((referenceItem, index) => (
              <MenuItem key={index} value={index}>{referenceItem}</MenuItem>
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
            value={category}
            onChange={handleCategoryChange}
            className={styles.select}
            select
            variant="standard"
          >
            {categories.map((categoryItem, index) => (
              <MenuItem key={index} value={index}>{categoryItem}</MenuItem>
            ))}
          </TextField>
          <Button variant="contained" className={styles.addButton} onClick={addCategory}>Add</Button>
        </Stack>
        <Stack
          direction="row"
          className={styles.tagStack}
        >
          {selectedIndexes.map((selectedIndex, index) => (
            <Tag key={index} text={categories[selectedIndex]} handleClick={() => removeTag(selectedIndex)} />
          ))}
        </Stack>
      </FormGroup>

    </Box>
  );
}
