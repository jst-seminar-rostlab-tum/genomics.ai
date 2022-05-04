import React, { useState } from 'react';
import {
  Stack, Button, Box, TextField, MenuItem, Typography,
} from '@mui/material';
import Tag from 'components/Tag';

// General filter that is needed in all categories
export default function FilterSelectTag(props) {
  const [category, setCategory] = useState('');
  const [selectedIndexes, setSelectedIndexes] = useState([]);
  const { categories, onChange, label } = props;

  const handleCategoryChange = (event) => {
    setCategory(event.target.value);
  };

  const handleIndexesAppend = (newValue) => {
    setSelectedIndexes((array) => [...array, newValue]);
    onChange(selectedIndexes);
  };

  const addCategory = () => {
    if (!selectedIndexes.includes(category) && category !== '') {
      handleIndexesAppend(category);
    }
  };

  const removeTag = (selectedIndex) => {
    setSelectedIndexes(selectedIndexes.filter((item) => item !== selectedIndex));
    onChange(selectedIndexes);
  };
  return (
    <Box>
      <Stack
        direction="row"
        spacing={2}
        alignItems="center"
        justifyContent="flex-end"
        sx={{ marginBottom: 2 }}
      >
        <Typography sx={{ fontWeight: 'bold' }}>{label}</Typography>
        <TextField
          value={category}
          onChange={handleCategoryChange}
          select
          variant="standard"
        >
          {categories.map((categoryItem, index) => (
            <MenuItem key={index} value={index}>{categoryItem}</MenuItem>
          ))}
        </TextField>
        <Button variant="contained" onClick={addCategory}>Add</Button>
      </Stack>
      <Stack
        direction="column"
        spacing={2}
      >
        {selectedIndexes.map((selectedIndex, index) => (
          <Tag
            key={index}
            content={categories[selectedIndex]}
            variant="primary-default"
            isDeletable
            handleDeletion={() => removeTag(selectedIndex)}
          />
        ))}
      </Stack>
    </Box>
  );
}
